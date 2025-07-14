# --- 1. Setup: Install and Load Required Packages ---
# This block ensures that all necessary packages are installed before use.
packages_to_install <- c("ggplot2", "dplyr", "tidyr", "xgboost", "viridis")
for (pkg in packages_to_install) {
  if (!require(pkg, character.only = TRUE)) {
    install.packages(pkg, dependencies = TRUE)
    library(pkg, character.only = TRUE)
  }
}

# --- 2. Create a Realistic Mock 'a_baumannii' DataFrame ---
# This section simulates the 'a_baumannii' dataframe you would have in your environment.
# It contains the essential columns: Year, Super_Region, and Cluster.
super_regions <- c(
  "East Asia and Pacific", "Eastern Europe and Central Asia", 
  "Latin America and Caribbean", "Middle East and North Africa", 
  "North America", "Eastern and Southern Africa", "Western Europe"
)

# Simulate individual isolate data
a_baumannii_mock <- expand_grid(
  Year = 2010:2023,
  Super_Region = super_regions,
  # Create a representative number of isolates per group
  isolate_id = 1:200 
) %>%
  mutate(
    # Assign clusters with varying probabilities
    random_prob = runif(n()),
    Cluster = case_when(
      random_prob < 0.4 + (Year - 2010) * 0.02 ~ "XDR-CRAB",
      random_prob < 0.6 ~ "MDR-CRAB",
      random_prob < 0.8 ~ "non-CRAB MDRAB",
      TRUE ~ "Pan-susceptible"
    )
  ) %>%
  select(-isolate_id, -random_prob)


# --- 3. Define the Master Plotting Function ---
# This single function takes the raw pathogen dataframe and encapsulates all 
# the logic for calculating proportions, forecasting, and plotting.

plot_xdr_crab_forecast <- function(a_baumannii_data, n_forecast_years = 5, n_bootstrap = 100) {
  
  # --- A. Calculate Proportions Internally ---
  # This makes the function self-contained.
  proportions_data <- a_baumannii_data %>%
    count(Super_Region, Year, Cluster) %>%
    group_by(Super_Region, Year) %>%
    mutate(proportion = n / sum(n)) %>%
    ungroup()
  
  # --- B. Nested Helper Function for Forecasting ---
  run_forecast_for_region <- function(region_data) {
    model_data <- region_data %>% arrange(Year)
    if(nrow(model_data) < 5) return(NULL)
    
    model_data_lagged <- model_data %>%
      mutate(proportion_lag1 = lag(proportion, 1)) %>%
      na.omit()
    
    if(nrow(model_data_lagged) < 2) return(NULL)
    
    initial_model <- xgboost(
      data = as.matrix(model_data_lagged$proportion_lag1),
      label = model_data_lagged$proportion,
      nrounds = 50, objective = "reg:squarederror", verbose = 0
    )
    residuals <- model_data_lagged$proportion - predict(initial_model, as.matrix(model_data_lagged$proportion_lag1))
    
    future_predictions_array <- array(NA, dim = c(n_bootstrap, n_forecast_years))
    
    for (i in 1:n_bootstrap) {
      bootstrap_target <- predict(initial_model, as.matrix(model_data_lagged$proportion_lag1)) + 
        sample(residuals, size = nrow(model_data_lagged), replace = TRUE)
      bootstrap_model <- xgboost(data = as.matrix(model_data_lagged$proportion_lag1), label = bootstrap_target, nrounds = 50, objective = "reg:squarederror", verbose = 0)
      
      last_known_value <- model_data %>% filter(Year == max(Year)) %>% pull(proportion)
      current_value <- last_known_value
      one_bootstrap_forecast <- numeric(n_forecast_years)
      for (j in 1:n_forecast_years) {
        prediction <- predict(bootstrap_model, matrix(current_value, nrow = 1))
        prediction <- max(0, min(1, prediction))
        one_bootstrap_forecast[j] <- prediction
        current_value <- prediction
      }
      future_predictions_array[i, ] <- one_bootstrap_forecast
    }
    
    forecast_df <- tibble(
      Year = (max(model_data$Year) + 1):(max(model_data$Year) + n_forecast_years),
      Point_Forecast = apply(future_predictions_array, 2, mean, na.rm = TRUE),
      LowerCI = apply(future_predictions_array, 2, quantile, probs = 0.025, na.rm = TRUE),
      UpperCI = apply(future_predictions_array, 2, quantile, probs = 0.975, na.rm = TRUE)
    )
    
    return(forecast_df)
  }
  
  # --- C. Main Function Logic ---
  
  crab_data <- proportions_data %>%
    filter(Cluster == "XDR-CRAB")
  
  regional_results <- crab_data %>%
    group_by(Super_Region) %>%
    nest() %>% 
    mutate(forecast = map(data, run_forecast_for_region)) %>%
    unnest(forecast) %>%
    select(Super_Region, Year, Point_Forecast, LowerCI, UpperCI)
  
  historicals_data <- crab_data %>%
    rename(Actual_Proportion = proportion)
  
  last_historical_points <- historicals_data %>%
    group_by(Super_Region) %>%
    filter(Year == max(Year)) %>%
    rename(Point_Forecast = Actual_Proportion)
  
  connected_forecast_data <- bind_rows(last_historical_points, regional_results) %>%
    arrange(Super_Region, Year)
  
  # --- D. Generate the final ggplot object ---
  final_plot <- ggplot() +
    geom_line(
      data = historicals_data,
      aes(x = Year, y = Actual_Proportion),
      color = "black", linewidth = 1
    ) +
    geom_ribbon(
      data = regional_results,
      aes(x = Year, ymin = LowerCI, ymax = UpperCI),
      fill = "skyblue", alpha = 0.5
    ) +
    geom_line(
      data = connected_forecast_data,
      aes(x = Year, y = Point_Forecast),
      color = "red", linetype = "dashed", linewidth = 1
    ) +
    facet_wrap(~ Super_Region, scales = "free_y") +
    scale_y_continuous(labels = scales::percent_format(accuracy = 1)) +
    scale_x_continuous(breaks = scales::pretty_breaks(n=4)) +
    labs(
      title = "Regional Forecast of XDR-CRAB Prevalence in A. baumannii",
      subtitle = "Historical Data and Forecast with 95% Confidence Interval",
      x = "Year",
      y = "Prevalence of XDR-CRAB Cluster"
    ) +
    theme_minimal(base_size = 14) +
    theme(
      strip.text = element_text(face = "bold"),
      axis.text.x = element_text(angle = 45, hjust = 1)
    )
  
  return(final_plot)
}


# --- 4. How to Use the Function ---

# Generate the plot by calling the function with your main data frame.
# In your real script, you would replace 'a_baumannii_mock' with 'a_baumannii'.
final_forecast_plot <- plot_xdr_crab_forecast(a_baumannii)

# Display the plot in RStudio
print(final_forecast_plot)

# (Optional) Save the plot to a file for your report
# ggsave("regional_xdr_crab_forecast.png", plot = final_forecast_plot, width = 14, height = 10, dpi = 300)
