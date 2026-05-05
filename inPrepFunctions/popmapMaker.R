library(DBI)
library(odbc)
library(dbplyr)
library(sf)
library(tidyverse)

# Define the path
barcode_path <- "C:/Users/jared/OneDrive - Michigan State University/MEER Lab - Documents/Projects/2025/2519_mueBaselineDevelopment/genotyping_2519/dependencies"

# List all .txt files
txt_files <- list.files(path = barcode_path, pattern = "barcodes\\.txt$", full.names = TRUE)

# Read and bind all files, keeping only the second column
sample_ids <- txt_files %>%
  map_dfr(~ read_tsv(.x, col_names = FALSE) %>% select(2)) %>%
  rename(sample_id = X2)



#Replace "C:/path/to/your/database.accdb" with the actual path
db_path <- "Z:/meerLab_database_be.accdb"

# Establish the connection
con <- dbConnect(odbc(),
                 Driver = "Microsoft Access Driver (*.mdb, *.accdb)",
                 Dbq = db_path)
# List tables
dbListTables(con)

samples <- tbl(con, "samples") %>% as_tibble()
collections <- tbl(con, "collections") %>% as_tibble()
locations <- tbl(con, "locations") %>% as_tibble()

files <- read_table(paste0(barcode_path, "/filenames.txt"),
                    col_names = FALSE)

left_join(samples, collections) %>% 
  left_join(locations) %>% 
  select(sample_id, location_name) %>% 
  filter(sample_id %in% sample_ids$sample_id) %>% 
  filter(sample_id %in% files$X1) %>% 
  write.table("C:/Users/jared/OneDrive - Michigan State University/MEER Lab - Documents/Projects/2025/2519_mueBaselineDevelopment/genotyping_2519/dependencies/pop.map",
              quote = FALSE,
              row.names = FALSE,
              col.names = FALSE)


# Close the connection
dbDisconnect(con)