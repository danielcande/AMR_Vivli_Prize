# --- 1. Load Required Packages ---

library(shiny)

library(tidyverse)

library(shinythemes)

library(sf)

library(rnaturalearth)
library(rnaturalearthdata)
library(countrycode)

library(xgboost)

library(gt) 

library(nnet)
data_env <- new.env()

load("vivli_data_compressed.RData", envir = data_env)

# --- 2. Check for and Combine Pre-Loaded Data ---

message("Starting app: Checking for required data objects...")

required_cluster_cols <- c(
  
  "a_baumannii" = "Cluster",
  
  "e_faecium"  = "Cluster",
  
  "e_spp"      = "Cluster",
  
  "s_aureus"   = "Cluster",
  
  "k_pneumoniae" = "Cluster",
  
  "p_aeruginosa" = "Cluster"
  
)


all_data_list <- lapply(names(required_cluster_cols), function(df_name) {
  
  if (!exists(df_name, envir = data_env)) {
    
    warning(paste("Data frame '", df_name, "' not found. Skipping."))
    
    return(NULL)
    
  }
  df <- get(df_name, envir = data_env) 
  
  actual_col_name <- "Cluster"
  
  if (actual_col_name %in% colnames(df)) {
    
    message(paste("Success: Found and processing", df_name))
    
    df <- df %>% dplyr::rename(Cluster = !!sym(actual_col_name))
    
    if (!"Species" %in% colnames(df)) {
      
      df$Species <- df_name
      
    }
    
    return(df)
    
  } else {
    
    warning(paste("Cluster column '", actual_col_name, "' not found in '", df_name, "'. Skipping."))
    
    return(NULL)
    
  }
  
})

all_data_processed <- bind_rows(all_data_list)

if(nrow(all_data_processed) == 0) {
  
  stop("App stopped: No data could be processed. Please check if your data frames are loaded and cluster column names are correct.")
  
}


# --- 3. Data Preparation ---

if (!"Super_Region" %in% colnames(all_data_processed)) all_data_processed$Super_Region <- "Unknown Region"

if (!"Country" %in% colnames(all_data_processed)) all_data_processed$Country <- "Unknown Country"

if (!"Isolate.ID" %in% colnames(all_data_processed)) all_data_processed$Isolate.ID <- 1:nrow(all_data_processed)

if (!"Year" %in% colnames(all_data_processed)) all_data_processed$Year <- 2022

if (!"Age.Group" %in% colnames(all_data_processed)) all_data_processed$Age.Group <- "Unknown"

if (!"Gender" %in% colnames(all_data_processed)) all_data_processed$Gender <- "Unknown"



world_map <- ne_countries(scale = "medium", returnclass = "sf") %>% dplyr::select(iso_a3, geometry)


all_data_processed <- all_data_processed %>% mutate(iso_a3 = countrycode(Country, "country.name", "iso3c"))

all_proportions_data <- all_data_processed %>%
  
  count(Species, Super_Region, Year, Cluster) %>%
  
  group_by(Species, Super_Region, Year) %>%
  
  mutate(proportion = n / sum(n)) %>%
  
  ungroup()

antibiotic_classes_lookup <- c(
  
  "Amikacin" = "Aminoglycoside", "Amoxycillin.clavulanate" = "Penicillin + Beta-lactamase inhibitor", "Ampicillin" = "Penicillin", "Cefepime" = "Cephalosporin (4th gen)", "Ceftazidime" = "Cephalosporin (3rd gen)",
  
  "Imipenem" = "Carbapenem", "Levofloxacin" = "Fluoroquinolone", "Meropenem" = "Carbapenem", "Piperacillin.tazobactam" = "Penicillin + Beta-lactamase inhibitor",
  
  "Tigecycline" = "Glycylcycline", "Ampicillin.sulbactam" = "Penicillin + Beta-lactamase inhibitor", "Aztreonam" = "Monobactam",
  
  "Cefixime" = "Cephalosporin (3rd gen)", "Ceftaroline" = "Cephalosporin (5th gen)", "Ceftazidime.avibactam" = "Cephalosporin + Beta-lactamase inhibitor",
  
  "Ciprofloxacin" = "Fluoroquinolone", "Colistin" = "Polymyxin", "Gentamicin" = "Aminoglycoside", "Trimethoprim.sulfa" = "Sulfonamide combination",
  
  "Ceftolozane.tazobactam" = "Cephalosporin + Beta-lactamase inhibitor", "Meropenem.vaborbactam" = "Carbapenem + Beta-lactamase inhibitor",
  
  "Cefpodoxime" = "Cephalosporin (3rd gen)", "Ceftibuten" = "Cephalosporin (3rd gen)",
  
  "Erythromycin" = "Macrolide", "Linezolid" = "Oxazolidinone", "Minocycline" = "Tetracycline",
  
  "Vancomycin" = "Glycopeptide", "Daptomycin" = "Lipopeptide", "Quinupristin.dalfopristin" = "Streptogramin",
  
  "Teicoplanin" = "Glycopeptide", "Oxacillin" = "Penicillin", "Clindamycin" = "Lincosamide",
  
  "Moxifloxacin" = "Fluoroquinolone", "Ceftriaxone" = "Cephalosporin (3rd gen)", "Doripenem" = "Carbapenem",
  
  "Ertapenem" = "Carbapenem"
  
)

antibiotic_classes_lookup <- antibiotic_classes_lookup[!duplicated(names(antibiotic_classes_lookup))]

antibiotic_classes_df <- tibble(Antibiotic = names(antibiotic_classes_lookup), Antibiotic_Class = unname(antibiotic_classes_lookup))

existing_antibiotics <- intersect(names(antibiotic_classes_lookup), colnames(all_data_processed))

all_data_final <- all_data_processed %>%
  
  pivot_longer(cols = all_of(existing_antibiotics), names_to = "Antibiotic", values_to = "Resistance") %>%
  
  left_join(antibiotic_classes_df, by = "Antibiotic") %>%
  
  filter(!is.na(Resistance), !is.na(Cluster))


all_data_final$Antibiotic_Class[is.na(all_data_final$Antibiotic_Class)] <- "Unknown Class"


age_choices <- c("All", unique(levels(all_data_final$Age.Group)))

gender_choices <- c("All", unique(levels(all_data_final$Gender)))



message("Data ready. Launching UI.")


# --- 4. Forecasting Function ---

generate_forecast_data <- function(proportions_data, target_species, target_region, n_forecast_years = 5, n_bootstrap = 100) {
  
  model_data <- proportions_data %>% filter(Species == target_species, Super_Region == target_region)
  
  feature_data <- model_data %>%
    
    dplyr::select(Year, Cluster, proportion) %>%
    
    pivot_wider(names_from = Cluster, values_from = proportion, names_prefix = "Cluster_", values_fill = 0) %>%
    
    arrange(Year)
  
  if(nrow(feature_data) < 10) {
    
    return(list(forecast_summary = NULL, historical_data = NULL, log = "Error: Insufficient historical data to generate a reliable forecast."))
    
  }
  
  feature_data_lagged <- feature_data %>% mutate(across(starts_with("Cluster_"), ~lag(.x, 1), .names = "{.col}_lag1"))
  
  predictor_variables <- colnames(feature_data_lagged)[grepl("_lag1$", colnames(feature_data_lagged))]
  
  target_variables <- colnames(feature_data)[grepl("^Cluster_", colnames(feature_data))]
  
  modeling_df <- feature_data_lagged %>% na.omit()
  
  all_bootstrap_models <- list()
  
  for (target_var in target_variables) {
    
    initial_model <- xgboost(data = as.matrix(modeling_df[, predictor_variables]), label = modeling_df[[target_var]], nrounds = 100, objective = "reg:squarederror", verbose = 0)
    
    residuals <- modeling_df[[target_var]] - predict(initial_model, as.matrix(modeling_df[, predictor_variables]))
    
    bootstrap_models_for_target <- list()
    
    for (i in 1:n_bootstrap) {
      
      bootstrap_target <- predict(initial_model, as.matrix(modeling_df[, predictor_variables])) + sample(residuals, size = nrow(modeling_df), replace = TRUE)
      
      bootstrap_model <- xgboost(data = as.matrix(modeling_df[, predictor_variables]), label = bootstrap_target, nrounds = 100, objective = "reg:squarederror", verbose = 0)
      
      bootstrap_models_for_target[[i]] <- bootstrap_model
      
    }
    
    all_bootstrap_models[[target_var]] <- bootstrap_models_for_target
    
  }
  
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
  
  forecast_summary_list <- list()
  
  for (k in 1:length(target_variables)) {
    
    target_var_name <- target_variables[k]
    
    predictions_for_one_cluster <- future_predictions_array[, , k]
    
    forecast_summary_list[[target_var_name]] <- data.frame(
      
      Year = (last_known_year + 1):(last_known_year + n_forecast_years),
      
      Point_Forecast = apply(predictions_for_one_cluster, 2, mean),
      
      LowerCI = apply(predictions_for_one_cluster, 2, quantile, probs = 0.025),
      
      UpperCI = apply(predictions_for_one_cluster, 2, quantile, probs = 0.975),
      
      Cluster = gsub("Cluster_", "", target_var_name)
      
    )
    
  }
  
  full_forecast_summary <- bind_rows(forecast_summary_list)
  
  plot_data_historical <- feature_data %>%
    
    pivot_longer(cols = starts_with("Cluster_"), names_to = "Cluster", values_to = "Actual_Proportion", names_prefix = "Cluster_")
  
  return(list(forecast_summary = full_forecast_summary, historical_data = plot_data_historical, log = "Forecast generation complete."))
  
}


# --- 5. Shiny UI (User Interface) ---

ui <- fluidPage(
  
  theme = shinytheme("cosmo"),
  
  titlePanel("AMR Pattern Explorer"),
  
  navbarPage("AMR Dashboard",
             
             tabPanel("Global Overview",
                      
                      fluidRow(column(12, wellPanel(h4("Shared Map Filters"), uiOutput("year_slider_global")))),
                      
                      fluidRow(
                        
                        column(6, wellPanel(h4("Global Cluster Map Filters"),
                                            
                                            selectInput("species_global_cluster", "Select Species:", choices = unique(all_data_final$Species)),
                                            
                                            selectInput("age_global_cluster", "Select Age Group:", choices = age_choices),
                                            
                                            radioButtons("gender_global_cluster", "Select Gender:", choices = gender_choices, inline = TRUE)
                                            
                        )),
                        
                        column(6, wellPanel(h4("Global Resistance Map Filters"),
                                            
                                            selectInput("species_global_res", "Select Species:", choices = unique(all_data_final$Species)),
                                            
                                            selectInput("analysis_level_global_res", "Analyse by:", choices = c("Antibiotic", "Antibiotic Class")),
                                            
                                            uiOutput("drug_select_global_res"),
                                            
                                            selectInput("age_global_res", "Select Age Group:", choices = age_choices),
                                            
                                            radioButtons("gender_global_res", "Select Gender:", choices = gender_choices, inline = TRUE)
                                            
                        ))
                        
                      ),
                      
                      fluidRow(column(6, plotOutput("cluster_map_global", height = "600px")), column(6, plotOutput("resistance_map_global", height = "600px")))
                      
             ),
             
             
             
             tabPanel("Resistance Profiles",
                      
                      sidebarLayout(
                        
                        sidebarPanel(
                          
                          h4("Profile Selection"),
                          
                          selectInput("species_heatmap", "Select Species to View Profile:",
                                      
                                      choices = c("A. baumannii" = "a_baumannii",
                                                  
                                                  "E. faecium" = "e_faecium",
                                                  
                                                  "E. spp" = "e_spp",
                                                  
                                                  "K. pneumoniae" = "k_pneumoniae",
                                                  
                                                  "P. aeruginosa" = "p_aeruginosa",
                                                  
                                                  "S. aureus" = "s_aureus"))
                          
                        ),
                        
                        mainPanel(
                          
                          gt_output(outputId = "resistance_heatmap_display")
                          
                        )
                        
                      )
                      
             ),
             
             tabPanel("Phenotype Resistance Trends",
                      
                      sidebarLayout(
                        
                        sidebarPanel(h4("Filter Options"), width = 3,
                                     
                                     selectInput("species_pheno", "1. Select Species:", choices = unique(all_data_final$Species)),
                                     
                                     selectInput("super_region_pheno", "2. Select Region:", choices = c("World", unique(levels(all_data_final$Super_Region)))),
                                     
                                     uiOutput("country_select_ui_pheno"), 
                                     
                                     selectInput("age_pheno", "4. Select Age Group:", choices = age_choices),
                                     
                                     radioButtons("gender_pheno", "5. Select Gender:", choices = gender_choices, inline = TRUE),
                                     
                                     uiOutput("year_slider_ui_pheno")),
                        
                        mainPanel(plotOutput("phenotype_plot", height = "600px"))
                        
                      )
                      
             ),
             
             tabPanel("Forecast Tool",
                      
                      sidebarLayout(
                        
                        sidebarPanel(h4("Forecast Model Controls"), width = 3, selectInput("species_forecast", "1. Select Species:", choices = unique(all_proportions_data$Species)), selectInput("region_forecast", "2. Select Super Region:", choices = unique(all_proportions_data$Super_Region)), actionButton("run_forecast_btn", "Run Forecast", icon = icon("chart-line")), hr(), uiOutput("cluster_select_ui_forecast")),
                        
                        mainPanel(plotOutput("forecast_plot", height = "600px"), h4("Model Log"), verbatimTextOutput("forecast_log"))
                        
                      )
                      
             ),
             
             tabPanel("Patient Risk Predictor",
                      
                      sidebarLayout(
                        
                        sidebarPanel(
                          
                          h4("Patient & Sample Details"),
                          
                          selectInput("pathogen_input", "Select Species of Suspected Infection:",
                                      
                                      choices = c("A. baumannii" = "a_baumannii",
                                                  
                                                  "E. faecium" = "e_faecium",
                                                  
                                                  "E. spp" = "e_spp",
                                                  
                                                  "K. pneumoniae" = "k_pneumoniae",
                                                  
                                                  "P. aeruginosa" = "p_aeruginosa",
                                                  
                                                  "S. aureus" = "s_aureus")),
                          
                          p("Enter the details below to predict the resistance profile."),
                          
                          
                          
                          radioButtons("pred_sexe", "Sex:", choices = c("Male", "Female"), inline = TRUE),
                          
                          
                          
                          selectInput("pred_dept", "Department:",
                                      
                                      choices = c("Nursing Home / Rehab", "Pediatric ICU", "Medicine General",
                                                  
                                                  "Surgery General", "Pediatric General", "Clinic / Office","Emergency Room",
                                                  
                                                  "Medicine ICU","None Given","Surgery ICU","General Unspecified","ICU Other")),
                          
                          
                          
                          selectInput("pred_region", "Region:",
                                      
                                      choices = levels(all_data_final$Super_Region)),
                          
                          selectInput("pred_age", "Age group:", choices = levels(all_data_final$Age.Group)),
                          
                          
                          
                          actionButton("predict_button", "Predict Profile", class = "btn-primary")
                          
                        ),
                        
                        mainPanel(
                          
                          h4("Predicted Probabilities"),
                          
                          p("The plot below shows the predicted probability of this isolate belonging to each resistance cluster."),
                          
                          plotOutput("prediction_plot")
                          
                        )
                        
                      )
                      
             )
             
  )
  
)


# --- 6. Shiny Server (Backend Logic) ---

server <- function(input, output, session) {
  
  # === Global Overview Tab ===
  
  output$year_slider_global <- renderUI({ year_range <- range(all_data_final$Year, na.rm = TRUE); sliderInput("year_range_global", "Select Year Range:", min = year_range[1], max = year_range[2], value = year_range, step = 1, sep = "") })
  
  output$drug_select_global_res <- renderUI({
    
    req(input$species_global_res, input$analysis_level_global_res)
    
    choices <- if (input$analysis_level_global_res == "Antibiotic") {
      
      all_data_final %>% filter(Species == input$species_global_res) %>% pull(Antibiotic) %>% unique() %>% sort()
      
    } else {
      
      all_data_final %>% filter(Species == input$species_global_res, !is.na(Antibiotic_Class)) %>% pull(Antibiotic_Class) %>% unique() %>% sort()
      
    }
    
    selectInput("drug_choice_global_res", "Select Antibiotic / Class:", choices = choices)
    
  })
  
  
  
  map_data_global_cluster <- reactive({
    
    req(input$species_global_cluster, input$year_range_global, input$age_global_cluster, input$gender_global_cluster)
    
    
    
    filtered_data <- all_data_final %>%
      
      filter(
        
        Species == input$species_global_cluster,
        
        Year >= input$year_range_global[1],
        
        Year <= input$year_range_global[2]
        
      )
    
    
    
    if(input$age_global_cluster != "All") {
      
      filtered_data <- filtered_data %>% filter(Age.Group == input$age_global_cluster)
      
    }
    
    if(input$gender_global_cluster != "All") {
      
      filtered_data <- filtered_data %>% filter(Gender == input$gender_global_cluster)
      
    }
    
    
    
    summary_data <- filtered_data %>%
      
      distinct(Isolate.ID, .keep_all = TRUE) %>%
      
      count(iso_a3, Cluster) %>%
      
      group_by(iso_a3) %>%
      
      slice_max(order_by = n, n = 1, with_ties = FALSE) %>%
      
      ungroup() %>%
      
      mutate(Dominant_Cluster = as.factor(Cluster))
    
    
    
    world_map %>% left_join(summary_data, by = "iso_a3")
    
  })
  
  
  
  map_data_global_resistance <- reactive({
    
    req(input$species_global_res, input$drug_choice_global_res, input$year_range_global, input$age_global_res, input$gender_global_res)
    
    
    
    filtered_data <- all_data_final %>%
      
      filter(
        
        Species == input$species_global_res,
        
        Year >= input$year_range_global[1],
        
        Year <= input$year_range_global[2],
        
        if (input$analysis_level_global_res == "Antibiotic") { Antibiotic == input$drug_choice_global_res } else { Antibiotic_Class == input$drug_choice_global_res }
        
      )
    
    
    
    if(input$age_global_res != "All") {
      
      filtered_data <- filtered_data %>% filter(Age.Group == input$age_global_res)
      
    }
    
    if(input$gender_global_res != "All") {
      
      filtered_data <- filtered_data %>% filter(Gender == input$gender_global_res)
      
    }
    
    
    
    summary_data <- filtered_data %>%
      
      group_by(iso_a3) %>%
      
      summarise(Prevalence = sum(Resistance == "Resistant", na.rm = TRUE) / n(), .groups = 'drop')
    
    
    
    world_map %>% left_join(summary_data, by = "iso_a3")
    
  })
  
  
  
  output$cluster_map_global <- renderPlot({ ggplot(data = map_data_global_cluster()) + geom_sf(aes(geometry = geometry, fill = Dominant_Cluster), color="white", linewidth=0.1) + scale_fill_viridis_d(na.value = "grey90", name = "Dominant Cluster") + labs(title = "Dominant Resistance Cluster by Country", subtitle = paste(input$species_global_cluster, "|", paste(input$year_range_global, collapse = "-"))) + theme_void() + theme(legend.position = "bottom") })
  
  output$resistance_map_global <- renderPlot({ ggplot(data = map_data_global_resistance()) + geom_sf(aes(geometry = geometry, fill = Prevalence), color="white", linewidth=0.1) + scale_fill_viridis_c(option = "magma", direction = -1, labels = scales::percent, na.value = "grey90", name = "Resistance") + labs(title = "Resistance Prevalence by Country", subtitle = paste(input$species_global_res, "|", input$drug_choice_global_res, "|", paste(input$year_range_global, collapse = "-"))) + theme_void() + theme(legend.position = "bottom") })
  
  
  
  output$resistance_heatmap_display <- render_gt({
    
    req(input$species_heatmap)
    
    switch(input$species_heatmap,
           
           "a_baumannii" = a_baumannii_heatmap,
           
           "e_faecium" = e_faecium_heatmap,
           
           "e_spp" = e_spp_heatmap,
           
           "k_pneumoniae" = k_pneumoniae_heatmap,
           
           "p_aeruginosa" = p_aeruginosa_heatmap,
           
           "s_aureus" = s_aureus_heatmap
           
    )
    
  })
  
  
  
  # === Phenotype Resistance Trends Tab ===
  
  
  
  output$country_select_ui_pheno <- renderUI({
    
    req(input$super_region_pheno)
    
    if (input$super_region_pheno == "World") {
      
      return(NULL)
      
    }
    
    # --- MODIFICATION: Explicitly convert Country to character to ensure correct values ---
    
    country_choices <- all_data_final %>%
      
      filter(Super_Region == input$super_region_pheno) %>%
      
      pull(Country) %>%
      
      as.character() %>% # This is the key fix
      
      unique() %>%
      
      sort()
    
    
    
    all_option <- paste("--- All of", input$super_region_pheno, "---")
    
    selectInput("country_pheno", "3. Select Country (Optional):", choices = c(all_option, country_choices))
    
  })
  
  
  
  output$year_slider_ui_pheno <- renderUI({
    
    year_range <- range(all_data_final$Year, na.rm = TRUE)
    
    sliderInput("year_range_pheno", "6. Select Year Range:", min = year_range[1], max = year_range[2], value = year_range, step = 1, sep = "")
    
  })
  
  
  
  filtered_data_pheno <- reactive({
    
    req(input$species_pheno, input$super_region_pheno, input$year_range_pheno, input$age_pheno, input$gender_pheno)
    
    
    
    base_data <- all_data_final %>%
      
      filter(
        
        Species == input$species_pheno,
        
        Year >= input$year_range_pheno[1],
        
        Year <= input$year_range_pheno[2]
        
      )
    
    
    
    if (input$age_pheno != "All") {
      
      base_data <- base_data %>% filter(Age.Group == input$age_pheno)
      
    }
    
    if (input$gender_pheno != "All") {
      
      base_data <- base_data %>% filter(Gender == input$gender_pheno)
      
    }
    
    
    
    if (input$super_region_pheno == "World") {
      
      return(base_data)
      
    } else {
      
      req(input$country_pheno)
      
      all_option_check <- paste("--- All of", input$super_region_pheno, "---")
      
      if (input$country_pheno == all_option_check) {
        
        base_data %>% filter(Super_Region == input$super_region_pheno)
        
      } else {
        
        base_data %>% filter(Country == input$country_pheno)
        
      }
      
    }
    
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
    
    
    
    # --- MODIFICATION: Corrected the location label logic ---
    
    location_label <- if (input$super_region_pheno == "World") {
      
      "World"
      
    } else {
      
      input$country_pheno # Directly use the input value
      
    }
    
    
    
    ggplot(plot_data_summary, aes(x = Year, y = Proportion, fill = as.factor(Cluster))) +
      
      geom_area(position = "fill", alpha = 0.8) +
      
      scale_y_continuous(labels = scales::percent_format()) +
      
      scale_x_continuous(breaks = scales::pretty_breaks()) +
      
      labs(
        
        title = "Temporal Trend of Resistance Clusters",
        
        subtitle = paste("Location:", location_label),
        
        x = "Year",
        
        y = "Proportion of Isolates",
        
        fill = "Cluster"
        
      ) +
      
      theme_minimal(base_size = 14) +
      
      theme(legend.position = "bottom")
    
  })
  
  
  
  # --- SERVER LOGIC FOR TAB 3: PREDICTION ---
  
  
  
  active_dataset <- reactive({
    
    req(input$pathogen_input)
    
    switch(input$pathogen_input,
           
           "a_baumannii" = a_baumannii,
           
           "e_faecium" = e_faecium,
           
           "e_spp" = e_spp,
           
           "p_aeruginosa" = p_aeruginosa,
           
           "s_aureus" = s_aureus,
           
           "k_pneumoniae" = k_pneumoniae)
    
  })
  
  
  
  prediction_result <- eventReactive(input$predict_button, {
    
    
    
    withProgress(message = 'Running Prediction', value = 0, {
      
      
      
      setProgress(value = 0.1, detail = "Loading data...")
      
      original_data <- active_dataset()
      
      
      
      predictor_formula <- Cluster ~ Super_Region + Gender + Age.Group + Speciality
      
      
      
      setProgress(value = 0.3, detail = "Training model...")
      
      model <- multinom(predictor_formula, data = original_data)
      
      
      
      new_data <- tibble(
        
        Age.Group = factor(input$pred_age, levels = levels(original_data$Age.Group)),
        
        Gender = factor(input$pred_sexe, levels = levels(original_data$Gender)),
        
        Speciality = factor(input$pred_dept, levels = levels(original_data$Speciality)),
        
        Super_Region = factor(input$pred_region, levels = levels(original_data$Super_Region))
        
      )
      
      
      
      setProgress(value = 0.6, detail = "Running bootstrap for CIs...")
      
      n_boot <- 20
      
      
      
      bootstrap_probs <- replicate(n_boot, {
        
        boot_sample <- original_data[sample(1:nrow(original_data), replace = TRUE), ]
        
        boot_model <- suppressWarnings(multinom(predictor_formula, data= boot_sample))
        
        predict(boot_model, newdata = new_data, type = "prob")
        
      })
      
      
      
      setProgress(value = 0.9, detail = "Finalizing results...")
      
      prob_cis <- apply(bootstrap_probs, 1, quantile, probs = c(0.025, 0.5, 0.975)) %>%
        
        t() %>%
        
        as.data.frame()
      
      
      
      colnames(prob_cis) <- c("LowerCI", "Probability", "UpperCI")
      
      prob_cis$Cluster <- rownames(prob_cis)
      
      
      
      return(prob_cis)
      
    }) # End withProgress
    
  })
  
  
  
  output$prediction_plot <- renderPlot({
    
    
    
    result_to_plot <- prediction_result()
    
    
    
    ggplot(result_to_plot, aes(x = Cluster, y = Probability, fill = Cluster)) +
      
      geom_col(alpha = 0.9) +
      
      geom_errorbar(aes(ymin = LowerCI, ymax = UpperCI), width = 0.2, linewidth = 0.8) +
      
      geom_text(aes(label = scales::percent(Probability, accuracy = 0.1)), vjust = -2.5, size = 5) +
      
      scale_y_continuous(limits = c(0, 1), labels = scales::percent) +
      
      scale_fill_brewer(palette = "Set2", guide = "none") +
      
      labs(
        
        x = "Predicted Resistance Cluster",
        
        y = "Probability (with 95% Confidence Interval)"
        
      ) +
      
      theme_minimal(base_size = 16)
    
  })
  
  
  
  # === Forecast Tool Tab ===
  
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
      
      selectInput("cluster_select_choice", "3. Select Cluster to Display:", choices = sort((choices)))
      
    }
    
  })
  
  output$forecast_plot <- renderPlot({
    
    req(forecast_data(), input$cluster_select_choice)
    
    results <- forecast_data()
    
    if (is.null(results$forecast_summary)) { return(ggplot() + labs(title = results$log) + theme_void()) }
    
    cluster_to_plot <- input$cluster_select_choice
    
    plot_data_hist <- results$historical_data %>% filter(Cluster == cluster_to_plot)
    
    plot_data_fcst <- results$forecast_summary %>% filter(Cluster == cluster_to_plot)
    
    if(nrow(plot_data_hist) == 0) return(ggplot() + labs(title=paste("No historical data for Cluster", cluster_to_plot)) + theme_void())
    
    last_actual_point <- plot_data_hist %>% filter(Year == max(Year))
    
    if(nrow(last_actual_point) > 0) {
      
      forecast_line_data <- bind_rows(data.frame(Year = last_actual_point$Year, Point_Forecast = last_actual_point$Actual_Proportion), plot_data_fcst %>% dplyr::select(Year, Point_Forecast))
      
    } else {
      
      forecast_line_data <- plot_data_fcst %>% dplyr::select(Year, Point_Forecast)
      
    }
    
    ggplot(plot_data_hist, aes(x = Year, y = Actual_Proportion)) +
      
      geom_line(aes(color = "Actual"), linewidth = 1.2) +
      
      geom_ribbon(data = plot_data_fcst, aes(x = Year, ymin = LowerCI, ymax = UpperCI), fill = "skyblue", alpha = 0.5, inherit.aes = FALSE) +
      
      geom_line(data = forecast_line_data, aes(x = Year, y = Point_Forecast, color = "Forecast"), linewidth = 1.2, linetype = "dashed") +
      
      labs(title = "XGBoost Forecast with 95% Confidence Interval", subtitle = paste("Forecast for Cluster", cluster_to_plot, "in", input$region_forecast), y = "Proportion of Isolates", color = "Legend") +
      
      scale_color_manual(values = c("Actual" = "black", "Forecast" = "red")) +
      
      theme_minimal(base_size = 14) +
      
      theme(legend.position = "bottom")
    
  })
  
  output$forecast_log <- renderText({ req(forecast_data()); forecast_data()$log })
  
}


# --- 7. Run the Shiny App ---

shinyApp(ui = ui, server = server) 
