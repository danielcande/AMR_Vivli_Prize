# --- 1. Load Required Packages ---
library(shiny)
library(tidyverse)
library(shinythemes)
library(sf)
library(rnaturalearth)
library(countrycode)
library(xgboost)
library(gt) 
library(nnet) 
library(DT)
# 
# #load("Vivli_final.RData")
# # At the top of app.R
load("/srv/shiny-server/AMR_Vivli_Prize/app_data.RData")

message("Data ready. Launching UI.")

# --- 4. Forecasting Function ---
generate_forecast_data <- function(proportions_data, target_species, target_region, n_forecast_years = 5, n_bootstrap = 20) {
  model_data <- proportions_data %>% dplyr::filter(Species == target_species, Super_Region == target_region)
  feature_data <- model_data %>%
    dplyr::select(Year, Cluster, proportion) %>%
    tidyr::pivot_wider(names_from = Cluster, values_from = proportion, names_prefix = "Cluster_", values_fill = 0) %>%
    dplyr::arrange(Year)
  if(nrow(feature_data) < 10) {
    return(list(forecast_summary = NULL, historical_data = NULL, log = "Error: Insufficient historical data to generate a reliable forecast."))
  }
  feature_data_lagged <- feature_data %>% dplyr::mutate(across(starts_with("Cluster_"), ~lag(.x, 1), .names = "{.col}_lag1"))
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
  last_known_features <- feature_data %>% dplyr::filter(Year == last_known_year) %>% dplyr::select(all_of(target_variables)) %>% as.numeric()
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
    tidyr::pivot_longer(cols = starts_with("Cluster_"), names_to = "Cluster", values_to = "Actual_Proportion", names_prefix = "Cluster_")
  return(list(forecast_summary = full_forecast_summary, historical_data = plot_data_historical, log = "Forecast generation complete."))
}

# --- 5. Enhanced Shiny UI with Medical Theme ---
ui <- fluidPage(
  # Custom CSS for medical/scientific theme
  tags$head(
    tags$style(HTML("
      /* Import medical-friendly fonts */
      @import url('https://fonts.googleapis.com/css2?family=Inter:wght@300;400;500;600;700&display=swap');
      
      /* Global body styling */
      body {
        font-family: 'Inter', 'Segoe UI', Tahoma, Geneva, Verdana, sans-serif;
        background: linear-gradient(135deg, #0f2027, #203a43, #2c5364); /* Darker, more professional gradient */
        min-height: 100vh;
      }
      
      /* --- FIX START --- */
      /* Make the main fluidPage container transparent */
      .container-fluid {
        background: transparent !important;
        box-shadow: none !important;
        padding-top: 0;
      }
      
      /* Style the CONTENT of the tabs as the white card */
      .tab-content > .tab-pane {
         background: rgba(255, 255, 255, 0.98);
         padding: 20px;
         margin-top: -20px; /* Pulls content up to meet the navbar */
         border-radius: 0 0 15px 15px;
         box-shadow: 0 10px 30px rgba(0, 0, 0, 0.2);
      }

      /* Make the main navbar transparent to show the body background */
      .navbar-default {
        background-color: transparent;
        border: none;
      }
      
      /* Change the default link color to white for legibility */
      .navbar-default .navbar-nav > li > a {
        color: white;
      }
      /* --- FIX END --- */
      
      /* Title styling */
      .navbar-brand {
        font-weight: 700 !important;
        font-size: 1.8rem !important;
        color: #2c3e50 !important;
        text-shadow: 0 2px 4px rgba(0,0,0,0.1);
      }
      
      /* Navigation tabs */
      .navbar-nav .nav-link {
        font-weight: 500;
        padding: 12px 20px !important;
        border-radius: 8px;
        margin: 0 3px;
        transition: all 0.3s ease;
      }
      
      .navbar-nav .nav-link:hover {
        background: linear-gradient(45deg, #1abc9c, #16a085); /* Teal hover */
        color: white !important;
        transform: translateY(-2px);
        box-shadow: 0 4px 12px rgba(26, 188, 156, 0.3);
      }
      
      .navbar-nav .nav-link.active, .navbar-default .navbar-nav > .active > a {
        background: linear-gradient(45deg, #3498db, #2980b9) !important; /* Blue active */
        color: white !important;
        box-shadow: 0 4px 15px rgba(52, 152, 219, 0.4);
      }
      
      /* Well panels - medical card style */
      .well {
        background: linear-gradient(145deg, #ffffff, #f8f9fa);
        border: 1px solid #e9ecef;
        border-radius: 12px;
        box-shadow: 0 4px 12px rgba(0, 0, 0, 0.08);
        margin-bottom: 20px;
        padding: 20px;
        transition: all 0.3s ease;
      }
      
      .well:hover {
        transform: translateY(-2px);
        box-shadow: 0 8px 25px rgba(0, 0, 0, 0.12);
      }
      
      /* Headers */
      h4 {
        color: #2c3e50;
        font-weight: 600;
        margin-bottom: 15px;
        border-bottom: 2px solid #3498db;
        padding-bottom: 8px;
      }
      
      /* Buttons */
      .btn-primary {
        background: linear-gradient(45deg, #3498db, #2980b9);
        border: none;
        border-radius: 8px;
        padding: 10px 20px;
        font-weight: 500;
        transition: all 0.3s ease;
        box-shadow: 0 4px 12px rgba(52, 152, 219, 0.3);
      }
      
      .btn-primary:hover {
        background: linear-gradient(45deg, #2980b9, #3498db);
        transform: translateY(-2px);
        box-shadow: 0 6px 20px rgba(52, 152, 219, 0.4);
      }
      
      .btn {
        border-radius: 8px;
        font-weight: 500;
        transition: all 0.3s ease;
      }
      
      /* Select inputs and controls */
      .form-control, .selectize-input {
        border-radius: 8px;
        border: 2px solid #e9ecef;
        transition: all 0.3s ease;
      }
      
      .form-control:focus, .selectize-input.focus {
        border-color: #3498db;
        box-shadow: 0 0 0 0.2rem rgba(52, 152, 219, 0.25);
      }
      
      /* Slider styling */
      .irs-bar {
        background: linear-gradient(45deg, #3498db, #2980b9);
      }
      
      .irs-handle {
        background: linear-gradient(45deg, #3498db, #2980b9);
        border-color: #2980b9;
      }
      
      /* Plot containers */
      .shiny-plot-output {
        border-radius: 12px;
        box-shadow: 0 4px 12px rgba(0, 0, 0, 0.1);
        overflow: hidden;
      }

      /* Guide and Glossary styling */
      .guide-section {
        background: linear-gradient(145deg, #ffffff, #f8f9fa);
        border-left: 4px solid #3498db;
        padding: 20px;
        margin: 15px 0;
        border-radius: 8px;
        box-shadow: 0 2px 8px rgba(0,0,0,0.1);
      }
      
      /* Verbatim output (logs) */
      pre {
        background: #2c3e50;
        color: #ecf0f1;
        border-radius: 8px;
        font-family: 'Consolas', 'Monaco', monospace;
        box-shadow: inset 0 2px 4px rgba(0,0,0,0.2);
      }

      /* --- NEW --- Disclaimer Box styling */
      .disclaimer-box {
          background-color: #f8d7da; /* Light red */
          color: #721c24; /* Dark red text */
          border: 1px solid #f5c6cb;
          border-radius: 8px;
          padding: 15px;
          margin: 15px 0;
          box-shadow: 0 2px 8px rgba(0,0,0,0.1);
      }
      .disclaimer-box h4 {
          color: #721c24 !important; /* Override default h4 color */
          border-bottom: none; /* Remove blue line */
          margin-top: 0;
      }
      /* --- END NEW --- */
    "))
  ),
  
  # Title Panel
  div(
    style = "text-align: center; padding: 20px 0; background: transparent; margin: 0; color: white;",
    h1("🧬 GUARDIAN", 
       style = "font-weight: 700; text-shadow: 0 2px 4px rgba(0,0,0,0.3); margin: 0;"),
    p("Global Understanding and Response Dashboard for Antimicrobial Resistance Intelligence and Navigation", 
      style = "font-size: 1.1rem; margin: 5px 0 0 0; font-weight: 300;")
  ),
  
  navbarPage(
    title = "", # Empty since we have custom title above
    theme = shinytheme("flatly"),
    # --- NEW WELCOME TAB ---
    # --- WELCOME TAB ---
    tabPanel("👋 Welcome",
             fluidPage(
               div(class = "well",
                   h2("Welcome to the GUARDIAN Dashboard", style = "text-align: center; color: #2c3e50;"),
                   br()
               ),
               
               div(class = "guide-section",
                   h4("About This Project"),
                   p("This dashboard provides an interactive platform for exploring global antimicrobial resistance (AMR) trends for the six ESKAPE pathogens. It is designed to help clinicians, public health officials, and researchers visualise complex surveillance data."),
                   p("This work was made possible by ", strong("Wellcome Trust"), " with access to critical data from ", strong("Pfizer"), " via the ", strong("Vivli"), " data sharing platform. We are grateful for their commitment to advancing AMR research.")
               ),
               
               # --- NEW SECTION ---
               div(class = "guide-section",
                   h4("🗺️ Navigating the Dashboard"),
                   p(strong("Resistance Profiles:"), " View detailed heatmaps that show the specific antibiotic susceptibility patterns for each identified resistance cluster."),
                   p(strong("Global Overview:"), " Explore interactive world maps showing the prevalence of different resistance clusters and specific antibiotic resistance rates."),
                   p(strong("Phenotype Resistance Trends:"), " Track how resistance to individual antibiotics has changed over time in different regions."),
                   p(strong("Forecast Tool:"), " Use a predictive model to forecast the prevalence of resistance clusters in the coming years."),
                   p(strong("Patient Risk Predictor:"), " Input patient details to predict the probability of their infection belonging to a specific resistance cluster."),
                   p(strong("How-To Guide:"), " Read case studies demonstrating how different users can leverage the dashboard's features."),
                   p(strong("Glossary:"), " A quick reference for the AMR-related acronyms and terms used throughout the dashboard.")
               ),
               # --- END NEW ---
               
               div(class = "guide-section",
                   h4("The ATLAS Dataset"),
                   p("The analyses are based on the Pfizer ", strong("ATLAS (Antimicrobial Testing Leadership and Surveillance)"), " database, a comprehensive global program collecting bacterial isolates. To prepare the data, we performed several cleaning steps:"),
                   tags$ul(
                     style = "margin-left: 20px;",
                     tags$li("Standardised antibiotic and species names."),
                     tags$li("Removed isolates with incomplete susceptibility profiles or critical missing metadata.")
                   )
               ),
               
               div(class = "guide-section",
                   h4("Important Considerations & Biases"),
                   p("Users should be aware of the inherent limitations of this dataset:"),
                   tags$ul(
                     style = "margin-left: 20px;",
                     tags$li(strong("Surveillance Bias:"), " Data is primarily from hospital settings and may not fully represent community-acquired infections."),
                     tags$li(strong("Geographic Bias:"), " Data collection is not uniform across all countries. Some regions are more densely sampled than others."),
                     tags$li(strong("Temporal Gaps:"), " Data availability may vary from year to year for specific locations.")
                   ),
                   p("This tool is intended for surveillance and research purposes only.")
               )
             )
    ),
    tabPanel("📖 How-To Guide",
             fluidPage(
               div(class = "well",
                   h2("📚 How to Use This Dashboard", style = "text-align: center; color: #2c3e50;"),
                   br()
               ),
               
               # --- NEW --- Disclaimer Box
               div(class = "disclaimer-box",
                   h4("⚠️ Disclaimer"),
                   p("This dashboard demonstrates a novel analytical approach using available surveillance data. For clinical decision support, implementation would require:"),
                   tags$ol(
                     style = "margin-left: 20px;", # Indent the list
                     tags$li("More representative local/regional datasets,"),
                     tags$li("Integration with clinical context including suspected infection aetiology, and"),
                     tags$li("Validation against patient outcomes.")
                   ),
                   p(strong("The methodology shown here provides a framework for such clinical applications."))
               ),
               # --- END NEW ---
               
               # --- Section for Clinicians ---
               div(class = "guide-section",
                   h4("🧑⚕️ For the Clinician at the Bedside"),
                   p("A 14-year-old boy arrives at an Emergency Room in New Delhi with a suspected soft tissue infection. The doctor needs to prescribe antibiotics now, days before lab results will be ready."),
                   p("The go-to choice is a broad-spectrum antibiotic, but the doctor first consults the ", strong("Patient Risk Predictor"), " tab. By entering the patient's details, the doctor sees a high probability that the infection belongs to an extensively drug-resistant (XDR) cluster."),
                   p(strong("Outcome:"), " The doctor is provided with a clinical support tool that can help to avoid a first-line treatment that would likely fail, instead exploring more targeted therapies and averting the risk of worsening the patient's condition.")
               ),
               
               # --- Section for Public Health Officers ---
               div(class = "guide-section",
                   h4("🔬 For the Public Health Surveillance Officer"),
                   p("An officer is alerted to a rumoured spike in carbapenem resistance in Northern England."),
                   p("They open the ", strong("Global Overview"), " tab and filter for the UK. They see that the dominant cluster in that region is indeed a known carbapenem-resistant type. Curious, they switch to the ", strong("Phenotype Resistance Trends"), " tab. Using the latest data, they confirm not only that carbapenem resistance is rising, but that this specific cluster has grown alongside it. A quick look at the ", strong("Forecast Tool"), " shows this trend is predicted to continue, highlighting the urgency of the issue."),
                   p(strong("Outcome:"), " The officer can quickly validate a public health concern, understand its scale, and use the forecast to inform regional infection control strategies and public health alerts.")
               ),
               
               # --- Section for Researchers ---
               div(class = "guide-section",
                   h4("👩🔬 For the Clinical Researcher"),
                   p("A researcher is planning their next grant proposal on Acinetobacter baumannii."),
                   p("They use the ", strong("Forecast Tool"), " to identify which resistance clusters are predicted to become most prevalent over the next five years. They notice a specific XDR cluster is expected to rise sharply. Intrigued, they examine its profile in the ", strong("Resistance Profiles"), " tab and see an unusual pattern of resistance to aminoglycosides. They go onto the ", strong("Patient Risk Predictor "),"tab and see how predicted membership in this cluster changes according to age, gender, location, and department."),
                   p(strong("Outcome:"), " The researcher has identified a high-impact, emerging resistance pattern to focus their research on, strengthening their grant application with data-driven evidence and ensuring their work addresses a future clinical need.")
               )
             )
    ),
    
    tabPanel("🧬 Glossary",
             fluidPage(
               div(class = "well",
                   h2("🔬 Glossary of AMR Terms", style = "text-align: center; color: #2c3e50;"),
                   br()
               ),
               
               # --- Resistance Categories ---
               div(class = "guide-section",
                   h4("🛡️ Resistance Categories"),
                   p(strong("MDR (Multi-Drug Resistant):"), " Bacteria that are resistant to at least one antibiotic in three or more different classes."),
                   p(strong("XDR (Extensively Drug-Resistant):"), " Bacteria that are resistant to all but two or fewer antibiotic classes, making them extremely difficult to treat.")
               ),
               
               # --- Specific Resistance Types ---
               div(class = "guide-section",
                   h4("🧪 Specific Resistance Types"),
                   p(strong("CRAB (Carbapenem-Resistant Acinetobacter baumannii):"), " A strain of Acinetobacter baumannii that is resistant to carbapenems, a class of last-resort antibiotics."),
                   p(strong("CRE (Carbapenem-Resistant Enterobacterales):"), " A family of bacteria (like E. coli and Klebsiella) that have developed resistance to carbapenem antibiotics."),
                   p(strong("CRP (Carbapenem-Resistant Pseudomonas):"), " A strain of Pseudomonas aeruginosa resistant to carbapenems, often found in hospital settings."),
                   p(strong("ESBL (Extended-Spectrum Beta-Lactamase):"), " Enzymes produced by some bacteria that allow them to break down and resist many common antibiotics, including penicillins and cephalosporins."),
                   p(strong("MRSA (Methicillin-Resistant Staphylococcus aureus):"), " A strain of Staphylococcus aureus that is resistant to methicillin and other related antibiotics, making it a common and challenging hospital-acquired infection.")
               ),
               
               # --- Infection Acquisition Types ---
               div(class = "guide-section",
                   h4("🏥 Infection Acquisition"),
                   p(strong("HA (Hospital-Acquired):"), " An infection that is contracted within a hospital or other healthcare facility."),
                   p(strong("CA (Community-Acquired):"), " An infection contracted by a person outside of a healthcare setting.")
               )
             )
    ),
    tabPanel("🔬 Resistance Profiles", 
             sidebarLayout( 
               sidebarPanel( 
                 div(class = "well",
                     h4("🧪 Profile Selection"), 
                     
                     # This selectInput now includes a 'selected' argument to ensure
                     # a heatmap loads by default when the tab is opened.
                     selectInput("species_heatmap", "Select Species to View Profile:", 
                                 choices = c("A. baumannii" = "a_baumannii", 
                                             "E. faecium" = "e_faecium", 
                                             "E. spp" = "e_spp", 
                                             "K. pneumoniae" = "k_pneumoniae", 
                                             "P. aeruginosa" = "p_aeruginosa", 
                                             "S. aureus" = "s_aureus"),
                                 selected = "a_baumannii"  # Default value
                     )
                 )
               ), 
               
               mainPanel( 
                 div(style = "padding: 20px;",
                     # Title and instructions for the new interactive table
                     h4("Resistance Prevalence by Class and Cluster"),
                     p("Click on any prevalence value to see a detailed breakdown by antibiotic.", style = "color: #555;"),
                     
                     # The DTOutput to display the interactive heatmap
                     DTOutput(outputId = "resistance_heatmap_display")
                 )
               ) 
             ) 
    ),
    
    tabPanel("🌍 Global Overview", 
             fluidRow(
               column(12, 
                      div(class = "well",
                          #style = "background: linear-gradient(135deg, #74b9ff, #0984e3);",
                          h4("🗺️ Shared Map Filters", style = "color: white; border-color: white;"), 
                          uiOutput("year_slider_global")
                      )
               )
             ), 
             fluidRow( 
               column(6, 
                      div(class = "well",
                          h4("🦠 Global Cluster Map Filters"), 
                          selectInput("species_global_cluster", "Select Species:", 
                                      choices = c("A. baumannii" = "a_baumannii", 
                                                  "E. faecium" = "e_faecium", 
                                                  "E. spp" = "e_spp", 
                                                  "K. pneumoniae" = "k_pneumoniae", 
                                                  "P. aeruginosa" = "p_aeruginosa", 
                                                  "S. aureus" = "s_aureus")), 
                          selectInput("age_global_cluster", "Select Age Group:", choices = age_choices), 
                          radioButtons("gender_global_cluster", "Select Gender:", choices = gender_choices, inline = TRUE) 
                      )
               ), 
               column(6, 
                      div(class = "well",
                          h4("💊 Global Resistance Map Filters"), 
                          selectInput("species_global_res", "Select Species:", 
                                      choices = c("A. baumannii" = "a_baumannii", 
                                                  "E. faecium" = "e_faecium", 
                                                  "E. spp" = "e_spp", 
                                                  "K. pneumoniae" = "k_pneumoniae", 
                                                  "P. aeruginosa" = "p_aeruginosa", 
                                                  "S. aureus" = "s_aureus")), 
                          selectInput("analysis_level_global_res", "Analyse by:", choices = c("Antibiotic", "Antibiotic Class")), 
                          uiOutput("drug_select_global_res"), 
                          selectInput("age_global_res", "Select Age Group:", choices = age_choices), 
                          radioButtons("gender_global_res", "Select Gender:", choices = gender_choices, inline = TRUE) 
                      )
               ) 
             ), 
             fluidRow(
               column(6, plotOutput("cluster_map_global", height = "600px")), 
               column(6, plotOutput("resistance_map_global", height = "600px"))
             ) 
    ), 
    
    
    tabPanel("📈 Phenotype Resistance Trends", 
             sidebarLayout( 
               sidebarPanel(
                 div(class = "well",
                     h4("⚙️ Filter Options"), 
                     selectInput("species_pheno", "1. Select Species:", 
                                 choices = c("A. baumannii" = "a_baumannii", 
                                             "E. faecium" = "e_faecium", 
                                             "E. spp" = "e_spp", 
                                             "K. pneumoniae" = "k_pneumoniae", 
                                             "P. aeruginosa" = "p_aeruginosa", 
                                             "S. aureus" = "s_aureus")), 
                     selectInput("super_region_pheno", "2. Select Region:", choices = c("World", unique(levels(all_data_final$Super_Region)))), 
                     uiOutput("country_select_ui_pheno"), 
                     selectInput("age_pheno", "4. Select Age Group:", choices = age_choices), 
                     radioButtons("gender_pheno", "5. Select Gender:", choices = gender_choices, inline = TRUE), 
                     uiOutput("year_slider_ui_pheno")
                 ),
                 width = 3
               ), 
               mainPanel(plotOutput("phenotype_plot", height = "600px")) 
             ) 
    ), 
    
    tabPanel("🔮 Forecast Tool", 
             sidebarLayout( 
               sidebarPanel(
                 div(class = "well",
                     h4("🤖 Forecast Model Controls"), 
                     selectInput("species_forecast", "1. Select Species:", 
                                 choices = c("A. baumannii" = "a_baumannii", 
                                             "E. faecium" = "e_faecium", 
                                             "E. spp" = "e_spp", 
                                             "K. pneumoniae" = "k_pneumoniae", 
                                             "P. aeruginosa" = "p_aeruginosa", 
                                             "S. aureus" = "s_aureus")), 
                     selectInput("region_forecast", "2. Select Super Region:", choices = unique(all_proportions_data$Super_Region)), 
                     actionButton("run_forecast_btn", "🚀 Run Forecast", icon = icon("chart-line"), class = "btn-primary"),
                     hr(), 
                     uiOutput("cluster_select_ui_forecast")
                 ),
                 width = 3
               ), 
               mainPanel(
                 plotOutput("forecast_plot", height = "600px"), 
                 div(class = "well",
                     h4("📋 Model Log"), 
                     verbatimTextOutput("forecast_log")
                 )
               ) 
             ) 
    ), 
    
    tabPanel("🏥 Patient Risk Predictor", 
             sidebarLayout( 
               sidebarPanel( 
                 div(class = "well",
                     h4("👤 Patient & Sample Details"), 
                     selectInput("pathogen_input", "Select Species of Suspected Infection:", 
                                 choices = c("A. baumannii" = "a_baumannii", 
                                             "E. faecium" = "e_faecium", 
                                             "E. spp" = "e_spp", 
                                             "K. pneumoniae" = "k_pneumoniae", 
                                             "P. aeruginosa" = "p_aeruginosa", 
                                             "S. aureus" = "s_aureus")), 
                     p("Enter the details below to predict the resistance profile.", style = "color: #7f8c8d;"), 
                     
                     radioButtons("pred_sexe", "👥 Gender:", choices = c("Male", "Female"), inline = TRUE), 
                     
                     selectInput("pred_dept", "🏢 Department:", 
                                 choices = c("Nursing Home / Rehab", "Pediatric ICU", "Medicine General", 
                                             "Surgery General", "Pediatric General", "Clinic / Office","Emergency Room", 
                                             "Medicine ICU","None Given","Surgery ICU","General Unspecified","ICU Other")), 
                     
                     selectInput("pred_region", "🌐 Region:", 
                                 choices = levels(all_data_final$Super_Region)), 
                     selectInput("pred_age", "📅 Age group:", choices = levels(k_pneumoniae$Age.Group)), 
                     
                     actionButton("predict_button", "🔍 Predict Profile", class = "btn-primary")
                 )
               ), 
               mainPanel( 
                 div(class = "well",
                     h4("📊 Predicted Probabilities"), 
                     p("The plot below shows the predicted probability of this isolate belonging to each resistance cluster.", 
                       style = "color: #7f8c8d;"), 
                     plotOutput("prediction_plot")
                 )
               ) 
             ) 
    )
  )
)

# --- 6. Shiny Server (Backend Logic) ---
server <- function(input, output, session) {
  # === Global Overview Tab ===
  output$year_slider_global <- renderUI({ 
    year_range <- range(all_data_final$Year, na.rm = TRUE)
    sliderInput("year_range_global", "Select Year Range:", min = year_range[1], max = year_range[2], value = year_range, step = 1, sep = "") 
  })
  
  output$drug_select_global_res <- renderUI({
    req(input$species_global_res, input$analysis_level_global_res)
    choices <- if (input$analysis_level_global_res == "Antibiotic") {
      all_data_final %>% dplyr::filter(Species == input$species_global_res) %>% pull(Antibiotic) %>% unique() %>% sort()
    } else {
      all_data_final %>% dplyr::filter(Species == input$species_global_res, !is.na(Antibiotic_Class)) %>% pull(Antibiotic_Class) %>% unique() %>% sort()
    }
    selectInput("drug_choice_global_res", "Select Antibiotic / Class:", choices = choices)
  })
  
  map_data_global_cluster <- reactive({
    req(input$species_global_cluster, input$year_range_global, input$age_global_cluster, input$gender_global_cluster)
    
    filtered_data <- all_data_final %>%
      dplyr::filter(
        Species == input$species_global_cluster,
        Year >= input$year_range_global[1],
        Year <= input$year_range_global[2]
      )
    
    if(input$age_global_cluster != "All") {
      filtered_data <- filtered_data %>% dplyr::filter(Age.Group == input$age_global_cluster)
    }
    if(input$gender_global_cluster != "All") {
      filtered_data <- filtered_data %>% dplyr::filter(Gender == input$gender_global_cluster)
    }
    
    # --- ADDITION: Calculate n ---
    n_isolates <- n_distinct(filtered_data$Isolate.ID)
    
    summary_data <- filtered_data %>%
      dplyr::distinct(Isolate.ID, .keep_all = TRUE) %>%
      dplyr::count(iso_a3, Cluster) %>%
      dplyr::group_by(iso_a3) %>%
      dplyr::slice_max(order_by = n, n = 1, with_ties = FALSE) %>%
      dplyr::ungroup() %>%
      dplyr::mutate(Dominant_Cluster = as.factor(Cluster))
    
    map_data <- world_map %>% dplyr::left_join(summary_data, by = "iso_a3")
    
    # --- MODIFICATION: Return a list with map data and the count ---
    return(list(map_data = map_data, n_isolates = n_isolates))
  })
  
  map_data_global_resistance <- reactive({
    # ... (similar filtering logic as above) ...
    req(input$species_global_res, input$drug_choice_global_res, input$year_range_global, input$age_global_res, input$gender_global_res)
    
    filtered_data <- all_data_final %>%
      dplyr::filter(
        Species == input$species_global_res,
        Year >= input$year_range_global[1],
        Year <= input$year_range_global[2],
        if (input$analysis_level_global_res == "Antibiotic") { Antibiotic == input$drug_choice_global_res } else { Antibiotic_Class == input$drug_choice_global_res }
      )
    
    if(input$age_global_res != "All") {
      filtered_data <- filtered_data %>% dplyr::filter(Age.Group == input$age_global_res)
    }
    if(input$gender_global_res != "All") {
      filtered_data <- filtered_data %>% dplyr::filter(Gender == input$gender_global_res)
    }
    
    # --- ADDITION: Calculate n ---
    n_isolates <- n_distinct(filtered_data$Isolate.ID)
    
    summary_data <- filtered_data %>%
      dplyr::group_by(iso_a3) %>%
      dplyr::summarise(
        Prevalence = sum(as.character(Resistance) == "Resistant", na.rm = TRUE) / n(), 
        .groups = 'drop'
      )
    
    map_data <- world_map %>% dplyr::left_join(summary_data, by = "iso_a3")
    
    # --- MODIFICATION: Return a list with map data and the count ---
    return(list(map_data = map_data, n_isolates = n_isolates))
  })
  
  output$cluster_map_global <- renderPlot({
    # --- MODIFICATION: Use the list and add caption ---
    plot_info <- map_data_global_cluster()
    n_caption <- paste("Total number of isolates:", format(plot_info$n_isolates, big.mark = ","))
    
    ggplot(data = plot_info$map_data) + 
      geom_sf(aes(geometry = geometry, fill = Dominant_Cluster), color="white", linewidth=0.1) + 
      scale_fill_viridis_d(na.value = "grey90", name = "Dominant Cluster") + 
      labs(title = "Dominant Resistance Cluster by Country", 
           subtitle = paste(input$species_global_cluster, "|", paste(input$year_range_global, collapse = "-")),
           caption = n_caption) + # Add caption
      theme_void() + 
      theme(legend.position = "bottom") 
  })
  
  output$resistance_map_global <- renderPlot({
    # --- MODIFICATION: Use the list and add caption ---
    plot_info <- map_data_global_resistance()
    n_caption <- paste("Total number of isolates:", format(plot_info$n_isolates, big.mark = ","))
    
    ggplot(data = plot_info$map_data) + 
      geom_sf(aes(geometry = geometry, fill = Prevalence), color="white", linewidth=0.1) + 
      scale_fill_viridis_c(option = "magma", direction = -1, labels = scales::percent, na.value = "grey90", name = "Resistance") + 
      labs(title = "Resistance Prevalence by Country", 
           subtitle = paste(input$species_global_res, "|", input$drug_choice_global_res, "|", paste(input$year_range_global, collapse = "-")),
           caption = n_caption) + # Add caption
      theme_void() + 
      theme(legend.position = "bottom") 
  })
  
  # 1. (NEW) Reactive expression to prepare data for the main heatmap table
  heatmap_data_reactive <- reactive({
    resistance_probs %>%
      dplyr::group_by(Antibiotic_Class, Class) %>%
      dplyr::summarise(Probability = mean(Probability, na.rm = TRUE), .groups = 'drop') %>%
      tidyr::pivot_wider(names_from = Class, values_from = Probability) %>%
      dplyr::arrange(Antibiotic_Class)
  })
  
  # 2. (MODIFIED) Replace render_gt with renderDT to create an interactive table
  heatmap_source_data <- reactive({
    req(input$species_heatmap)
    
    # 1. First, select the appropriate pre-loaded 'gt_tbl' object based on user input
    gt_object <- switch(input$species_heatmap,
                        "a_baumannii" = a_baumannii_heatmap,
                        "e_faecium"   = e_faecium_heatmap,
                        "e_spp"       = e_spp_heatmap,
                        "k_pneumoniae"  = k_pneumoniae_heatmap,
                        "p_aeruginosa"  = p_aeruginosa_heatmap,
                        "s_aureus"    = s_aureus_heatmap
    )
    
    # 2. Extract the underlying data frame from the selected gt object.
    #    This is the key step that fixes the error.
    data_frame <- gt_object$`_data`
    
    return(data_frame)
  })
  
  # 2. (MODIFIED) Render the interactive DT table from the selected source data.
  output$resistance_heatmap_display <- renderDT({
    
    dt_data <- heatmap_source_data()
    cluster_cols <- colnames(dt_data)[-1] # Get cluster columns for formatting
    
    datatable(
      dt_data,
      class = 'cell-border stripe hover',
      selection = 'none',
      rownames = FALSE,
      options = list(
        dom = 't', 
        ordering = FALSE, 
        pageLength = nrow(dt_data),
        scrollX = TRUE),
      # This JavaScript sends cell click info back to Shiny
      callback = JS(
        "table.on('click', 'td:not(:first-child)', function() {",
        "  var class_name = $(this).closest('tr').find('td:first').text();",
        "  var col_idx = $(this).index();",
        "  var cluster_name = $(this).closest('table').find('th').eq(col_idx).text();",
        "  if(class_name && cluster_name) {",
        "    Shiny.setInputValue('heatmap_cell_clicked', { class: class_name, cluster: cluster_name, nonce: Math.random() });",
        "  }",
        "});"
      )
    ) %>%
      formatStyle(
        columns = cluster_cols,
        backgroundColor = styleInterval(c(0.2, 0.4, 0.6, 0.8), c('#FEF0D9', '#FDCC8A', '#FC8D59', '#E34A33', '#B30000')),
        color = styleInterval(0.5, c('black', 'white'))
      ) %>%
      formatPercentage(columns = cluster_cols, digits = 0)
  })
  
  # 3. (UNCHANGED) The observer for the drill-down modal.
  # This part still uses your global 'resistance_probs' object as requested.
  observeEvent(input$heatmap_cell_clicked, {
    
    click_info <- input$heatmap_cell_clicked
    
    # --- FIX: Swap the variables to match your data structure ---
    # The value from the column header (the Antibiotic Class) is in 'click_info$cluster'
    # The value from the row name (the Resistance Cluster) is in 'click_info$class'
    class_name <- click_info$cluster
    cluster_name <- click_info$class
    
    # --- (You can remove the debugging print statements now) ---
    
    drilldown_data <- resistance_probs %>%
      dplyr::filter(
        Antibiotic_Class == class_name,
        Class == cluster_name
      )
    
    # (The rest of your code is correct and does not need to be changed)
    drilldown_plot <- ggplot(drilldown_data, aes(x = factor(Antibiotic, levels = Antibiotic), y = Probability, fill = Probability)) +
      geom_col(show.legend = FALSE) +
      geom_text(aes(label = scales::percent(Probability, accuracy = 1)), hjust = -0.2, size = 4) +
      coord_flip() +
      scale_fill_viridis_c(option = "magma", direction = -1) +
      scale_y_continuous(labels = scales::percent, limits = c(0, 1.1), breaks = seq(0, 1, 0.2)) +
      labs(
        title = paste("Resistance within Class:", class_name),
        subtitle = paste("Showing all species for Cluster:", cluster_name),
        x = NULL,
        y = "Prevalence of Resistance"
      ) +
      theme_minimal(base_size = 14) +
      theme(panel.grid.major.y = element_blank(), panel.grid.minor.x = element_blank(), plot.title = element_text(face = "bold"))
    
    showModal(modalDialog(
      title = "Detailed Resistance Profile",
      renderPlot(drilldown_plot),
      easyClose = TRUE,
      footer = modalButton("Close"),
      size = "l"
    ))
  })
  
  # === Phenotype Resistance Trends Tab ===
  output$country_select_ui_pheno <- renderUI({
    req(input$super_region_pheno)
    if (input$super_region_pheno == "World") {
      return(NULL)
    }
    
    country_choices <- all_data_final %>% 
      dplyr::filter(Super_Region == input$super_region_pheno) %>% 
      pull(Country) %>% 
      unique() %>% 
      as.character() %>%
      sort()  
    
    all_option <- paste("--- All of", input$super_region_pheno, "---")
    selectInput("country_pheno", "3. Select Country (Optional):", 
                choices = c(all_option, country_choices))
  })
  
  output$year_slider_ui_pheno <- renderUI({
    year_range <- range(all_data_final$Year, na.rm = TRUE)
    sliderInput("year_range_pheno", "6. Select Year Range:", min = year_range[1], max = year_range[2], value = year_range, step = 1, sep = "")
  })
  
  filtered_data_pheno <- reactive({
    req(input$species_pheno, input$super_region_pheno, input$year_range_pheno, input$age_pheno, input$gender_pheno)
    
    base_data <- all_data_final %>%
      dplyr::filter(
        Species == input$species_pheno,
        Year >= input$year_range_pheno[1],
        Year <= input$year_range_pheno[2]
      )
    
    if (input$age_pheno != "All") {
      base_data <- base_data %>% dplyr::filter(Age.Group == input$age_pheno)
    }
    if (input$gender_pheno != "All") {
      base_data <- base_data %>% dplyr::filter(Gender == input$gender_pheno)
    }
    
    if (input$super_region_pheno == "World") {
      return(base_data)
    } else {
      req(input$country_pheno)
      all_option_check <- paste("--- All of", input$super_region_pheno, "---")
      if (input$country_pheno == all_option_check) {
        base_data %>% dplyr::filter(Super_Region == input$super_region_pheno)
      } else {
        base_data %>% dplyr::filter(Country == input$country_pheno)
      }
    }
  })
  
  output$phenotype_plot <- renderPlot({
    plot_data <- filtered_data_pheno()
    if (nrow(plot_data) == 0) return(ggplot() + labs(title = "No data available for this selection") + theme_void())
    
    # --- ADDITION: Calculate n and create caption ---
    n_isolates <- n_distinct(plot_data$Isolate.ID)
    n_caption <- paste("Total number of isolates:", format(n_isolates, big.mark = ","))
    
    plot_data_summary <- plot_data %>%
      dplyr::distinct(Isolate.ID, .keep_all = TRUE) %>%
      dplyr::count(Year, Cluster) %>%
      dplyr::group_by(Year) %>%
      dplyr::mutate(Proportion = n / sum(n)) %>%
      dplyr::ungroup()
    
    location_label <- if (input$super_region_pheno == "World") {
      "World"
    } else {
      input$country_pheno
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
        fill = "Cluster",
        caption = n_caption # Add caption
      ) +
      theme_minimal(base_size = 14) +
      theme(legend.position = "bottom")
  })
  
  # === Patient Risk Predictor Tab ===
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
      model <- multinom(predictor_formula, data = original_data, trace = FALSE)
      
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
        boot_model <- suppressWarnings(multinom(predictor_formula, data = boot_sample, trace = FALSE))
        predict(boot_model, newdata = new_data, type = "prob")
      })
      
      setProgress(value = 0.9, detail = "Finalizing results...")
      prob_cis <- apply(bootstrap_probs, 1, quantile, probs = c(0.025, 0.5, 0.975)) %>%
        t() %>%
        as.data.frame()
      
      colnames(prob_cis) <- c("LowerCI", "Probability", "UpperCI")
      prob_cis$Cluster <- rownames(prob_cis)
      
      return(prob_cis)
    })
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
    
    plot_data_hist <- results$historical_data %>% dplyr::filter(Cluster == cluster_to_plot)
    
    plot_data_fcst <- results$forecast_summary %>% dplyr::filter(Cluster == cluster_to_plot)
    
    if(nrow(plot_data_hist) == 0) return(ggplot() + labs(title=paste("No historical data for Cluster", cluster_to_plot)) + theme_void())
    
    last_actual_point <- plot_data_hist %>% dplyr::filter(Year == max(Year))
    
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
