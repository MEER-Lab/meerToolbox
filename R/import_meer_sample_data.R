#' Fetch and Align Metadata for GTscore Samples
#' 
#' @param gen_obj The filtered genind object
#' @param db_path Path to the MS Access database (Default: MEER Lab Backend)
#' @param select_cols A character vector of column names. Defaults to everything().
#'
import_meer_sample_data <- function(gen_obj, 
                                 db_path = "Z:/meerLab_database_be.accdb", 
                                 select_cols = NULL) {
  
  # 1. Load Dependencies
  if (!requireNamespace("odbc", quietly = TRUE)) stop("Package 'odbc' required.")
  if (!requireNamespace("DBI", quietly = TRUE)) stop("Package 'DBI' required.")
  library(tidyverse)
  library(DBI)
  
  # 2. Database Connection
  con <- dbConnect(odbc::odbc(),
                   Driver = "Microsoft Access Driver (*.mdb, *.accdb)",
                   Dbq = db_path)
  
  # 3. Import Tables
  samples <- tbl(con, "samples") %>% as_tibble()
  collections <- tbl(con, "collections") %>% as_tibble()
  locations <- tbl(con, "locations") %>% as_tibble()
  dbDisconnect(con)
  
  # 4. Join Logic
  message("Joining tables and aligning with genind samples...")
  expected_names <- adegenet::indNames(gen_obj)
  
  # Join using suffixes to prevent column collisions
  metadata_full <- samples %>%
    filter(sample_id %in% expected_names) %>%
    left_join(collections, by = "collection_id", suffix = c("", "_coll")) %>%
    left_join(locations, by = "location_id", suffix = c("", "_loc"))
  
  # 5. Strict Sample Count Validation
  if (nrow(metadata_full) != length(expected_names)) {
    missing <- setdiff(expected_names, metadata_full$sample_id)
    stop(paste0("CRITICAL ERROR: Sample count mismatch!\n",
                "Genind has ", length(expected_names), " samples, but DB returned ", nrow(metadata_full), ".\n",
                "Check these samples in the database: ", paste(missing, collapse = ", ")))
  }
  
  # 6. Alignment and Selection
  # Default behavior: if select_cols is NULL, use everything()
  target_cols <- if(is.null(select_cols)) colnames(metadata_full) else select_cols
  
  # Safety check: only select columns that actually exist
  available_cols <- intersect(target_cols, colnames(metadata_full))
  
  metadata_final <- metadata_full %>%
    arrange(match(sample_id, expected_names)) %>%
    select(all_of(available_cols))
  
  # 7. Final Order Validation
  if (!all(metadata_final$sample_id == expected_names)) {
    stop("CRITICAL ERROR: Sample data alignment failed. Sample order does not match genind.")
  }
  
  message(paste("Success: Returned", nrow(metadata_final), "samples and", ncol(metadata_final), "columns."))
  return(as_tibble(metadata_final))
}