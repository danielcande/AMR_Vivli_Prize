# --- create_database.R ---

# Step 1: Install necessary packages
# install.packages(c("DBI", "RSQLite", "dplyr", "tidyr"))

library(DBI)
library(RSQLite)
library(dplyr)
library(tidyr)

# Step 2: Load your existing .RData file
load("vivli_data.RData") 

# Step 3: Create a connection to a new SQLite database
con <- dbConnect(RSQLite::SQLite(), "vivli_database.sqlite")

# Step 4: Get a list of all the data frames to save
data_frame_names <- c(
  "a_baumannii", "e_faecium", "e_spp", 
  "s_aureus", "k_pneumoniae", "p_aeruginosa"
  # Add any other core data frames here
)

# Step 5: Combine all data frames into one large table
# This is more efficient for querying later
all_data_list <- lapply(data_frame_names, function(df_name) {
  df <- get(df_name)
  if (!"Species" %in% colnames(df)) {
    df$Species <- df_name
  }
  return(df)
})
all_data_combined <- bind_rows(all_data_list)

# Write the combined data to a single table
dbWriteTable(con, "all_data_combined", all_data_combined, overwrite = TRUE)
cat("Wrote main table: all_data_combined\n")

# Step 6: Create and store the proportions data needed for the forecast tab
all_proportions_data <- all_data_combined %>%
  count(Species, Super_Region, Year, Cluster) %>%
  group_by(Species, Super_Region, Year) %>%
  mutate(proportion = n / sum(n)) %>%
  ungroup()

dbWriteTable(con, "all_proportions_data", all_proportions_data, overwrite = TRUE)
cat("Wrote summary table: all_proportions_data\n")

# Step 7: Store data for the heatmaps
# Instead of storing complex gt objects, we store the data needed to build them.
for (df_name in data_frame_names) {
  df <- get(df_name)
  # Assuming your heatmap data comes from a summary like this:
  heatmap_data <- df %>%
    group_by(Cluster) %>%
    summarise(across(where(is.numeric), ~mean(.x, na.rm = TRUE))) # Example summary
  
  dbWriteTable(con, paste0(df_name, "_heatmap_data"), heatmap_data, overwrite = TRUE)
  cat("Wrote heatmap data for:", df_name, "\n")
}


# Step 8: Close the connection
dbDisconnect(con)

cat("\nDatabase 'vivli_database.sqlite' created successfully!\n")

