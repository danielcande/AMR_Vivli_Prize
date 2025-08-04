# --- 1. Load Required Packages ---
# You need these libraries to run the processing steps.
library(tidyverse)
library(sf)
library(rnaturalearth)
library(countrycode)

# --- 2. Load Raw Data ---
# Load the original, large .RData file that contains all the raw data frames.
message("Loading raw data from Vivli_final.RData...")
load("Vivli_final.RData") # Make sure this file is in your working directory

# --- 3. Perform All Slow Data Processing ---
message("Starting heavy data processing...")

# Combine the individual data frames
required_cluster_cols <- c(
  "a_baumannii" = "Cluster",
  "e_faecium"   = "Cluster",
  "e_spp"       = "Cluster",
  "s_aureus"    = "Cluster",
  "k_pneumoniae" = "Cluster",
  "p_aeruginosa" = "Cluster"
)

# Clean up any bad data before combining
e_faecium <- e_faecium %>%
  dplyr::filter(Gender != "-")

all_data_list <- lapply(names(required_cluster_cols), function(df_name) {
  if (!exists(df_name, envir = .GlobalEnv)) {
    warning(paste("Data frame '", df_name, "' not found. Skipping."))
    return(NULL)
  }
  df <- get(df_name, envir = .GlobalEnv)
  
  if ("Cluster" %in% colnames(df)) {
    message(paste("Success: Found and processing", df_name))
    if (!"Species" %in% colnames(df)) {
      df$Species <- df_name
    }
    return(df)
  } else {
    warning(paste("Cluster column not found in '", df_name, "'. Skipping."))
    return(NULL)
  }
})

all_data_processed <- bind_rows(all_data_list)

# Add missing columns if they don't exist
if (!"Super_Region" %in% colnames(all_data_processed)) all_data_processed$Super_Region <- "Unknown Region"
if (!"Country" %in% colnames(all_data_processed)) all_data_processed$Country <- "Unknown Country"
if (!"Isolate.ID" %in% colnames(all_data_processed)) all_data_processed$Isolate.ID <- 1:nrow(all_data_processed)
if (!"Year" %in% colnames(all_data_processed)) all_data_processed$Year <- 2022
if (!"Age.Group" %in% colnames(all_data_processed)) all_data_processed$Age.Group <- "Unknown"
if (!"Gender" %in% colnames(all_data_processed)) all_data_processed$Gender <- "Unknown"

# Create the final, long-format data frame for the app
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
  dplyr::mutate(iso_a3 = countrycode(Country, "country.name", "iso3c")) %>%
  tidyr::pivot_longer(cols = all_of(existing_antibiotics), names_to = "Antibiotic", values_to = "Resistance") %>%
  dplyr::left_join(antibiotic_classes_df, by = "Antibiotic") %>%
  dplyr::filter(!is.na(Resistance), !is.na(Cluster))

all_data_final$Antibiotic_Class[is.na(all_data_final$Antibiotic_Class)] <- "Unknown Class"

# Create other necessary objects for the app
all_proportions_data <- all_data_processed %>%
  dplyr::count(Species, Super_Region, Year, Cluster) %>%
  dplyr::group_by(Species, Super_Region, Year) %>%
  dplyr::mutate(proportion = n / sum(n)) %>%
  dplyr::ungroup()

world_map <- ne_countries(scale = "medium", returnclass = "sf") %>% 
  dplyr::select(iso_a3, geometry)

age_choices <- c("All", unique(all_data_final$Age.Group))
gender_choices <- c("All", unique(all_data_final$Gender))


# --- 4. Save Final Objects ---
# Save only the clean, final objects that the Shiny app needs to run.
message("Saving final objects to app_data.RData...")
save(
  # The main data frames needed by the app
  all_data_final,
  all_proportions_data,
  world_map,
  age_choices,
  gender_choices,
  
  # The individual data frames needed for the predictor tab
  a_baumannii,
  e_faecium,
  e_spp,
  s_aureus,
  k_pneumoniae,
  p_aeruginosa,
  
  # The heatmap objects needed for the profiles tab
  # NOTE: This assumes these heatmap objects exist in your Vivli_final.RData file.
  # If they are created differently, you'll need to add that code to this script.
  a_baumannii_heatmap,
  e_faecium_heatmap,
  e_spp_heatmap,
  k_pneumoniae_heatmap,
  p_aeruginosa_heatmap,
  s_aureus_heatmap,
  
  # The file to save to
  file = "app_data.RData"
)

message("Done! 'app_data.RData' is ready to be uploaded to the server.")

