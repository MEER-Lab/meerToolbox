
library(DBI)
library(odbc)
library(dbplyr)
library(sf)
library(tidyverse)

#Replace "C:/path/to/your/database.accdb" with the actual path
db_path <- "C:/Users/jared/OneDrive/Desktop/MeerLabDB.accdb"

# Establish the connection
con <- dbConnect(odbc(),
                 Driver = "Microsoft Access Driver (*.mdb, *.accdb)",
                 Dbq = db_path)
# List tables
dbListTables(con)

samples <- tbl(con, "samples") %>% as_tibble()
collections <- tbl(con, "collections") %>% as_tibble()
locations <- tbl(con, "locations") %>% as_tibble()

# View the icollections# View the imported data
head(data_from_access_dbi)

# Close the connection
dbDisconnect(con)