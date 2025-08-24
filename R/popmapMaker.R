library(DBI)
library(odbc)
library(dbplyr)
library(sf)
library(tidyverse)



# Define the path
barcode_path <- "C:/Users/jared/OneDrive - Michigan State University/Shared Documents - MEER Lab/Projects/2025/2507_KsSdWaeSbfI/genotyping_2507_2521/dependencies"

# List all .txt files
txt_files <- list.files(path = barcode_path, pattern = "\\.txt$", full.names = TRUE)

# Read and bind all files, keeping only the second column
sample_ids <- txt_files %>%
  map_dfr(~ read_tsv(.x, col_names = FALSE) %>% select(2)) %>%
  rename(sample_id = X2)



#Replace "C:/path/to/your/database.accdb" with the actual path
db_path <- "C:/Users/jared/OneDrive/Desktop/meerLab_database.accdb"

# Establish the connection
con <- dbConnect(odbc(),
                 Driver = "Microsoft Access Driver (*.mdb, *.accdb)",
                 Dbq = db_path)
# List tables
dbListTables(con)

samples <- tbl(con, "samples") %>% as_tibble()
collections <- tbl(con, "collections") %>% as_tibble()
locations <- tbl(con, "locations") %>% as_tibble()

left_join(samples, collections) %>% 
  left_join(locations) %>% 
  select(sample_id, location_name) %>% 
  filter(sample_id %in% sample_ids$sample_id)


# Close the connection
dbDisconnect(con)