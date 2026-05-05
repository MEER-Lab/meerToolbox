library(readxl)
library(tidyverse)

# Path to your Excel file
excel_file <- "C:/Users/jared/OneDrive - Michigan State University/MEER Lab - Documents/Projects/2025/2523_waeBroodstock20242025/2523 & 2509 Walleye gt-seq barcodes v2.xlsx"

# Get all sheet names
sheets <- excel_sheets(excel_file)

# Read the first sheet:
# - Skip the first row
# - Use the second row as header
first_sheet <- read_excel(
  path = excel_file,
  sheet = sheets[1],
  skip = 1) %>%
  mutate(sheet_name = sheets[1])

# Read the remaining sheets:
# - No rows to skip
# - No headers, so set col_names = FALSE
other_sheets <- map(sheets[-1], ~ {
  read_excel(
    path = excel_file,
    sheet = .x,
    col_names = FALSE) %>%
    mutate(PlateID = .x)}) %>%
  bind_rows() %>% 
  rename("sampleBarcode" = 1,
         "sampleID" = 2)

# Combine everything into one tibble
left_join(first_sheet, other_sheets) %>% 
  mutate(barcode = paste0(i7, "-", sampleBarcode)) %>% 
  select(sampleID, barcode) %>% 
  write_delim("C:/Users/jared/OneDrive - Michigan State University/MEER Lab - Documents/Projects/2025/2523_waeBroodstock20242025/genotyping_2523/barcodes25092523.txt",
              delim = "\t",
              col_names = FALSE)
