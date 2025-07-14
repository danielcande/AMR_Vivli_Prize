# --- create_database.R ---

# Step 1: Install the necessary packages if you haven't already
# install.packages(c("DBI", "RSQLite"))

library(DBI)
library(RSQLite)
library(dplyr)

# Step 2: Load your existing .RData file into the environment
# This will load all your data frames (a_baumannii, e_faecium, etc.)
load("vivli_data.RData") 

# Step 3: Create a connection to a new SQLite database file
# This will create the file "vivli_database.sqlite" in your project folder.
con <- dbConnect(RSQLite::SQLite(), "vivli_database.sqlite")

# Step 4: Get a list of all the data frame objects you just loaded
# We assume they are the main data frames in your environment.
# Be sure to add any other data frames you need to this list.
data_frame_names <- c(
  "a_baumannii", 
  "e_faecium", 
  "e_spp", 
  "s_aureus", 
  "k_pneumoniae", 
  "p_aeruginosa"
  # Add any other data frames that are in your .RData file here
)

# Step 5: Write each data frame to a table in the database
for (df_name in data_frame_names) {
  # Get the actual data frame object from its name
  df_object <- get(df_name)
  
  # Write the data frame to the database. 
  # The table name will be the same as the data frame name.
  dbWriteTable(con, df_name, df_object, overwrite = TRUE)
  
  # Print a progress message
  cat("Wrote table:", df_name, "\n")
}

# Step 6: Close the connection
dbDisconnect(con)

cat("\nDatabase 'vivli_database.sqlite' created successfully!\n")
