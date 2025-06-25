# --- AMR Resistance Profile Dashboard with Integrated Forecast Tool ---
# This Shiny app visualizes AMR data and includes a tab for
# running a predictive XGBoost model to forecast future cluster trends.

# --- 1. Load Required Packages ---
library(shiny)
library(tidyverse)
library(shinythemes)
library(sf)
library(rnaturalearth)
library(countrycode)
library(xgboost)

# --- 2. Check for and Combine Pre-Loaded Data ---
message("Starting app: Checking for required data objects...")

required_cluster_cols <- c(
  "a_baumannii" = "a_baumannii_model$predclass",
  "e_faecium"   = "e_faecium_model$predclass",
  "e_spp"       = "e_spp_model$predclass",
  "s_aureus"    = "s_aureus_model$predclass"
)

all_data_list <- lapply(names(required_cluster_cols), function(df_name) {
  if (!exists(df_name, envir = .GlobalEnv)) {
    warning(paste("Data frame '", df_name, "' not found. Skipping."))
    return(NULL)
  }
  df <- get(df_name, envir = .GlobalEnv)
  cluster_col_name <- required_cluster_cols[[df_name]]
  if (cluster_col_name %in% colnames(df)) {
    message(paste("Success: Found and processing", df_name))
    df <- df %>% dplyr::rename(Cluster = !!sym(cluster_col_name))
    if (!"Species" %in% colnames(df)) {
      df$Species <- df_name
    }
    return(df)
  } else {
    warning(paste("Cluster column '", cluster_col_name, "' not found in '", df_name, "'. Skipping."))
    return(NULL)
  }
})

all_data_processed <- bind_rows(all_data_list)

if(nrow(all_data_processed) == 0) {
  stop("App stopped: No data could be processed. Please check if your data frames are loaded and cluster column names are correct.")
}

# --- Add default columns for robustness ---
if (!"Super_Region" %in% colnames(all_data_processed)) all_data_processed$Super_Region <- "Unknown Region"
if (!"Country" %in% colnames(all_data_processed)) all_data_processed$Country <- "Unknown Country"
if (!"Isolate.ID" %in% colnames(all_data_processed)) all_data_processed$Isolate.ID <- 1:nrow(all_data_processed)
if (!"Year" %in% colnames(all_data_processed)) all_data_processed$Year <- 2022

# --- Prepare Geospatial Data ---
world_map <- ne_countries(scale = "medium", returnclass = "sf") %>%
  dplyr::select(iso_a3, geometry) # Use standard ISO code

# Standardize country names in your data
all_data_processed <- all_data_processed %>%
  mutate(iso_a3 = countrycode(Country, "country.name", "iso3c"))

# --- Prepare Proportional Data for Forecasting ---
all_proportions_data <- all_data_processed %>%
  count(Species, Super_Region, Year, Cluster) %>%
  group_by(Species, Super_Region, Year) %>%
  mutate(proportion = n / sum(n)) %>%
  ungroup()

# --- Prepare Reshaped Data for Other Plots ---
antibiotic_classes_lookup <- c(
  "Amikacin" = "Aminoglycoside", "Amoxycillin.clavulanate" = "Penicillin + Beta-lactamase inhibitor", "Ampicillin" = "Penicillin", "Cefepime" = "Cephalosporin (4th gen)", "Ceftazidime" = "Cephalosporin (3rd gen)", "Imipenem" = "Carbapenem", "Levofloxacin" = "Fluoroquinolone", "Meropenem" = "Carbapenem", "Piperacillin.tazobactam" = "Penicillin + Beta-lactamase inhibitor", "Tigecycline" = "Glycylcycline", "Ampicillin.sulbactam" = "Penicillin + Beta-lactamase inhibitor", "Aztreonam" = "Monobactam", "Cefixime" = "Cephalosporin (3rd gen)", "Ceftaroline" = "Cephalosporin (5th gen)", "Ceftazidime.avibactam" = "Cephalosporin + Beta-lactamase inhibitor", "Ciprofloxacin" = "Fluoroquinolone", "Colistin" = "Polymyxin", "Gentamicin" = "Aminoglycoside", "Trimethoprim.sulfa" = "Sulfonamide combination", "Ceftolozane.tazobactam" = "Cephalosporin + Beta-lactamase inhibitor", "Meropenem.vaborbactam" = "Carbapenem + Beta-lactamase inhibitor", "Cefpodoxime" = "Cephalosporin (3rd gen)", "Ceftibuten" = "Cephalosporin (3rd gen)", "Erythromycin" = "Macrolide", "Linezolid" = "Oxazolidinone", "Minocycline" = "Tetracycline", "Vancomycin" = "Glycopeptide", "Daptomycin" = "Lipopeptide", "Quinupristin.dalfopristin" = "Streptogramin", "Teicoplanin" = "Glycopeptide", "Oxacillin" = "Penicillin", "Clindamycin" = "Lincosamide", "Moxifloxacin" = "Fluoroquinolone"
)
antibiotic_classes_lookup <- antibiotic_classes_lookup[!duplicated(names(antibiotic_classes_lookup))] # Remove duplicates

antibiotic_classes_df <- tibble(Antibiotic = names(antibiotic_classes_lookup), Antibiotic_Class = unname(antibiotic_classes_lookup))
existing_antibiotics <- intersect(names(antibiotic_classes_lookup), colnames(all_data_processed))
all_data_final <- all_data_processed %>%
  pivot_longer(cols = all_of(existing_antibiotics), names_to = "Antibiotic", values_to = "Resistance") %>%
  left_join(antibiotic_classes_df, by = "Antibiotic") %>%
  filter(!is.na(Resistance), !is.na(Cluster))
all_data_final$Antibiotic_Class[is.na(all_data_final$Antibiotic_Class)] <- "Unknown Class"

message("Data ready. Launching UI.")

# ===================================================================
# --- Forecasting Function ---
# ===================================================================
generate_forecast_data <- function(proportions_data, target_species, target_region, n_forecast_years = 5, n_bootstrap = 50) {
  
  # 1. Feature Engineering
  model_data <- proportions_data %>%
    filter(Species == target_species, Super_Region == target_region)
  
  feature_data <- model_data %>%
    dplyr::select(Year, Cluster, proportion) %>%
    pivot_wider(names_from = Cluster, values_from = proportion, names_prefix = "Cluster_", values_fill = 0) %>%
    arrange(Year)
  
  if(nrow(feature_data) < 10) {
    return(list(forecast_summary = NULL, historical_data = NULL, log = "Error: Insufficient historical data to generate a reliable forecast."))
  }
  
  feature_data_lagged <- feature_data %>%
    mutate(across(starts_with("Cluster_"), ~lag(.x, 1), .names = "{.col}_lag1"))
  
  # 2. Multi-Model Training with Bootstrapping
  predictor_variables <- colnames(feature_data_lagged)[grepl("_lag1$", colnames(feature_data_lagged))]
  target_variables <- colnames(feature_data)[grepl("^Cluster_", colnames(feature_data))]
  modeling_df <- feature_data_lagged %>% na.omit()
  all_bootstrap_models <- list()
  
  for (target_var in target_variables) {
    initial_model <- xgboost(data = as.matrix(modeling_df[, predictor_variables]), label = modeling_df[[target_var]], nrounds = 50, objective = "reg:squarederror", verbose = 0)
    residuals <- modeling_df[[target_var]] - predict(initial_model, as.matrix(modeling_df[, predictor_variables]))
    bootstrap_models_for_target <- list()
    for (i in 1:n_bootstrap) {
      bootstrap_target <- predict(initial_model, as.matrix(modeling_df[, predictor_variables])) + sample(residuals, size = nrow(modeling_df), replace = TRUE)
      bootstrap_model <- xgboost(data = as.matrix(modeling_df[, predictor_variables]), label = bootstrap_target, nrounds = 50, objective = "reg:squarederror", verbose = 0)
      bootstrap_models_for_target[[i]] <- bootstrap_model
    }
    all_bootstrap_models[[target_var]] <- bootstrap_models_for_target
  }
  
  # 3. Iterative Forecasting
  last_known_year <- max(feature_data$Year)
  last_known_features <- feature_data %>% filter(Year == last_known_year) %>% dplyr::select(all_of(target_variables)) %>% as.numeric()
  future_predictions_array <- array(NA, dim = c(n_bootstrap, n_forecast_years, length(target_variables)))
  
  for (i in 1:n_bootstrap) {
    current_features <- last_known_features
    for (j in 1:n_forecast_years) {
      yearly_predictions <- numeric(length(target_variables))
      input_matrix <- matrix(current_features, nrow = 1); colnames(input_matrix) <- predictor_variables
      for (k in 1:length(target_variables)) {
        target_var <- target_variables[k]
        current_model <- all_bootstrap_models[[target_var]][[i]]
        yearly_predictions[k] <- predict(current_model, input_matrix)
      }
      yearly_predictions[yearly_predictions < 0] <- 0
      normalized_predictions <- yearly_predictions / sum(yearly_predictions, na.rm=T)
      future_predictions_array[i, j, ] <- normalized_predictions
      current_features <- normalized_predictions
    }
  }
  
  # 4. Summarize Forecasts
  forecast_summary_list <- list()
  for (k in 1:length(target_variables)) {
    target_var_name <- target_variables[k]
    predictions_for_one_cluster <- future_predictions_array[, , k]
    forecast_summary_list[[target_var_name]] <- data.frame(
      Year = (last_known_year + 1):(last_known_year + n_forecast_years),
      Point_Forecast = apply(predictions_for_one_cluster, 2, mean),
      Lower_CI = apply(predictions_for_one_cluster, 2, quantile, probs = 0.025),
      Upper_CI = apply(predictions_for_one_cluster, 2, quantile, probs = 0.975),
      Cluster = gsub("Cluster_", "", target_var_name)
    )
  }
  
  full_forecast_summary <- bind_rows(forecast_summary_list)
  plot_data_historical <- feature_data %>%
    pivot_longer(cols = starts_with("Cluster_"), names_to = "Cluster", values_to = "Actual_Proportion", names_prefix = "Cluster_")
  
  return(list(forecast_summary = full_forecast_summary, historical_data = plot_data_historical, log = "Forecast generation complete."))
}


# --- 3. Shiny UI (User Interface) ---
ui <- fluidPage(
  theme = shinytheme("cosmo"),
  titlePanel("AMR Pattern Explorer"),
  navbarPage("AMR Dashboard",
             tabPanel("Global Overview",
                      fluidRow(
                        column(6, wellPanel(
                          h4("Global Cluster Map Filters"),
                          selectInput("species_global_cluster", "Select Species:", choices = unique(all_data_final$Species)),
                          uiOutput("year_slider_global_cluster")
                        )),
                        column(6, wellPanel(
                          h4("Global Resistance Map Filters"),
                          selectInput("species_global_res", "Select Species:", choices = unique(all_data_final$Species)),
                          selectInput("analysis_level_global_res", "Analyze by:", choices = c("Antibiotic", "Antibiotic Class")),
                          uiOutput("drug_select_global_res")
                        ))
                      ),
                      fluidRow(
                        column(6, plotOutput("cluster_map_global", height = "600px")),
                        column(6, plotOutput("resistance_map_global", height = "600px"))
                      )
             ),
             tabPanel("Phenotype Resistance Trends",
                      sidebarLayout(
                        sidebarPanel(
                          h4("Filter Options"), width = 3,
                          selectInput("species_pheno", "1. Select Species:", choices = unique(all_data_final$Species)),
                          selectInput("super_region_pheno", "2. Select Super Region:", choices = unique(all_data_final$Super_Region)),
                          uiOutput("country_select_ui_pheno"),
                          uiOutput("year_slider_ui_pheno")
                        ),
                        mainPanel(plotOutput("phenotype_plot", height = "600px"))
                      )
             ),
             tabPanel("Forecast Tool",
                      sidebarLayout(
                        sidebarPanel(
                          h4("Forecast Model Controls"), width = 3,
                          selectInput("species_forecast", "1. Select Species:", choices = unique(all_proportions_data$Species)),
                          selectInput("region_forecast", "2. Select Super Region:", choices = unique(all_proportions_data$Super_Region)),
                          actionButton("run_forecast_btn", "Run Forecast", icon = icon("chart-line")),
                          hr(),
                          # This UI will be populated after the model runs
                          uiOutput("cluster_select_ui_forecast")
                        ),
                        mainPanel(
                          plotOutput("forecast_plot", height = "600px"),
                          h4("Model Log"),
                          verbatimTextOutput("forecast_log")
                        )
                      ))
  )
)

# --- 4. Shiny Server (Backend Logic) ---
server <- function(input, output, session) {
  
  # ==================================================
  # === Global Overview Tab: Server Logic ====
  # ==================================================
  output$year_slider_global_cluster <- renderUI({ year_range <- range(all_data_final$Year, na.rm = TRUE); sliderInput("year_range_global_cluster", "Select Year Range:", min = year_range[1], max = year_range[2], value = year_range, step = 1, sep = "") })
  output$drug_select_global_res <- renderUI({ req(input$species_global_res, input$analysis_level_global_res); choices <- if (input$analysis_level_global_res == "Antibiotic") { all_data_final %>% filter(Species == input$species_global_res) %>% pull(Antibiotic) %>% unique() %>% sort() } else { all_data_final %>% filter(Species == input$species_global_res, !is.na(Antibiotic_Class)) %>% pull(Antibiotic_Class) %>% unique() %>% sort() }; selectInput("drug_choice_global_res", "Select Antibiotic / Class:", choices = choices) })
  
  map_data_global_cluster <- reactive({
    req(input$species_global_cluster, input$year_range_global_cluster)
    summary_data <- all_data_final %>% 
      filter(Species == input$species_global_cluster, Year >= input$year_range_global_cluster[1], Year <= input$year_range_global_cluster[2]) %>% 
      distinct(Isolate.ID, .keep_all = TRUE) %>% 
      count(iso_a3, Cluster) %>% 
      group_by(iso_a3) %>% 
      slice_max(order_by = n, n = 1, with_ties = FALSE) %>% 
      ungroup() %>% 
      mutate(Dominant_Cluster = as.factor(Cluster))
    world_map %>% left_join(summary_data, by = "iso_a3")
  })
  
  map_data_global_resistance <- reactive({
    req(input$species_global_res, input$drug_choice_global_res)
    summary_data <- all_data_final %>% 
      filter(Species == input$species_global_res, 
             if (input$analysis_level_global_res == "Antibiotic") { 
               Antibiotic == input$drug_choice_global_res 
             } else { 
               Antibiotic_Class == input$drug_choice_global_res 
             }) %>% 
      group_by(iso_a3) %>% 
      summarise(Prevalence = sum(Resistance == "Resistant", na.rm = TRUE) / n(), .groups = 'drop')
    world_map %>% left_join(summary_data, by = "iso_a3")
  })
  
  output$cluster_map_global <- renderPlot({
    ggplot(data = map_data_global_cluster()) +
      geom_sf(aes(geometry = geometry, fill = Dominant_Cluster), color="white", size=0.1) +
      scale_fill_viridis_d(na.value = "grey90", name = "Dominant Cluster") +
      labs(title = "Dominant Resistance Cluster by Country", subtitle = paste(input$species_global_cluster, "|", paste(input$year_range_global_cluster, collapse = "-"))) +
      theme_void() + theme(legend.position = "bottom")
  })
  
  output$resistance_map_global <- renderPlot({
    ggplot(data = map_data_global_resistance()) +
      geom_sf(aes(geometry = geometry, fill = Prevalence), color="white", size=0.1) +
      scale_fill_viridis_c(option = "magma", direction = -1, labels = scales::percent, na.value = "grey90", name = "Resistance") +
      labs(title = "Resistance Prevalence by Country", subtitle = paste(input$species_global_res, "|", input$drug_choice_global_res)) +
      theme_void() + theme(legend.position = "bottom")
  })
  
  # ==================================================
  # === Phenotype Resistance Trends Tab: Server Logic (UPDATED) ===
  # ==================================================
  output$country_select_ui_pheno <- renderUI({ req(input$super_region_pheno); country_choices <- all_data_final %>% filter(Super_Region == input$super_region_pheno) %>% pull(Country) %>% unique() %>% sort(); all_option <- paste("--- All of", input$super_region_pheno, "---"); selectInput("country_pheno", "3. Select Country (Optional):", choices = c(all_option, country_choices)) })
  output$year_slider_ui_pheno <- renderUI({ year_range <- range(all_data_final$Year, na.rm = TRUE); sliderInput("year_range_pheno", "4. Select Year Range:", min = year_range[1], max = year_range[2], value = year_range, step = 1, sep = "") })
  
  filtered_data_pheno <- reactive({
    req(input$species_pheno, input$country_pheno, input$year_range_pheno)
    all_option_check <- paste("--- All of", input$super_region_pheno, "---")
    
    geo_filtered <- if (input$country_pheno == all_option_check) {
      all_data_final %>% filter(Super_Region == input$super_region_pheno)
    } else {
      all_data_final %>% filter(Country == input$country_pheno)
    }
    
    geo_filtered %>%
      filter(Species == input$species_pheno, Year >= input$year_range_pheno[1], Year <= input$year_range_pheno[2])
  })
  
  output$phenotype_plot <- renderPlot({
    plot_data <- filtered_data_pheno()
    if (nrow(plot_data) == 0) return(ggplot() + labs(title = "No data available for this selection") + theme_void())
    
    plot_data_summary <- plot_data %>%
      distinct(Isolate.ID, .keep_all = TRUE) %>%
      count(Year, Cluster) %>%
      group_by(Year) %>%
      mutate(Proportion = n / sum(n)) %>%
      ungroup()
    
    ggplot(plot_data_summary, aes(x = Year, y = Proportion, fill = as.factor(Cluster))) +
      geom_area(position = "fill", alpha = 0.8) +
      scale_y_continuous(labels = scales::percent_format()) +
      scale_x_continuous(breaks = scales::pretty_breaks()) +
      labs(title = "Temporal Trend of Resistance Clusters", subtitle = paste("Location:", input$country_pheno),
           x = "Year", y = "Proportion of Isolates", fill = "Cluster") +
      theme_minimal(base_size = 14) +
      theme(legend.position = "bottom")
  })
  
  # ============================================
  # === Forecast Tool Tab: Server Logic ========
  # ============================================
  
  forecast_data <- eventReactive(input$run_forecast_btn, {
    req(input$species_forecast, input$region_forecast)
    withProgress(message = 'Training models and generating forecast...', value = 0, {
      results <- generate_forecast_data(proportions_data = all_proportions_data, target_species = input$species_forecast, target_region = input$region_forecast)
      return(results)
    })
  })
  
  output$cluster_select_ui_forecast <- renderUI({
    req(forecast_data())
    results <- forecast_data()
    if(!is.null(results$forecast_summary)){
      choices <- unique(results$forecast_summary$Cluster)
      selectInput("cluster_forecast_choice", "3. Select Cluster to Display:", choices = sort((choices)))
    }
  })
  
  output$forecast_plot <- renderPlot({
    req(forecast_data(), input$cluster_forecast_choice)
    results <- forecast_data()
    
    if (is.null(results$forecast_summary)) {
      return(results$plot) # Show error plot from function
    }
    
    cluster_to_plot <- input$cluster_forecast_choice
    
    plot_data_hist <- results$historical_data %>% filter(Cluster == cluster_to_plot)
    plot_data_fcst <- results$forecast_summary %>% filter(Cluster == cluster_to_plot)
    
    if(nrow(plot_data_hist) == 0) return(ggplot() + labs(title=paste("No historical data for Cluster", cluster_to_plot)) + theme_void())
    
    last_actual_point <- plot_data_hist %>% filter(Year == max(Year))
    
    if(nrow(last_actual_point) > 0) {
      forecast_line_data <- bind_rows(
        data.frame(Year = last_actual_point$Year, Point_Forecast = last_actual_point$Actual_Proportion),
        plot_data_fcst %>% dplyr::select(Year, Point_Forecast)
      )
    } else {
      forecast_line_data <- plot_data_fcst %>% dplyr::select(Year, Point_Forecast)
    }
    
    ggplot(plot_data_hist, aes(x = Year, y = Actual_Proportion)) +
      geom_line(aes(color = "Actual"), linewidth = 1.2) +
      geom_ribbon(data = plot_data_fcst, aes(x = Year, ymin = Lower_CI, ymax = Upper_CI), fill = "skyblue", alpha = 0.5, inherit.aes = FALSE) +
      geom_line(data = forecast_line_data, aes(x = Year, y = Point_Forecast, color = "Forecast"), linewidth = 1.2, linetype = "dashed") +
      labs(title = "XGBoost Forecast with 95% Confidence Interval", subtitle = paste("Forecast for Cluster", cluster_to_plot, "in", input$region_forecast), y = "Proportion of Isolates", color = "Legend") +
      scale_color_manual(values = c("Actual" = "black", "Forecast" = "red")) +
      theme_minimal(base_size = 14) +
      theme(legend.position = "bottom")
  })
  
  output$forecast_log <- renderText({
    req(forecast_data())
    forecast_data()$log
  })
  
}

# --- 5. Run the Shiny App ---
shinyApp(ui = ui, server = server)
