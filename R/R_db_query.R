
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

collections %>% 
  filter(species_common == "cisco") %>% 
  select(species_common, location_name, total_samples) %>% 
  group_by(location_name) %>% 
  summarize(sum_samples = sum(total_samples)) %>% 
  write.table("C:/Users/jared/OneDrive/Desktop/ciscoSamplesAtMSU.txt",
              quote = FALSE,
              row.names = FALSE)
  
spatDat <- collections %>% 
  filter(species_common == "cisco") %>%
  group_by(location_name) %>% 
  summarize(sum_samples = sum(total_samples)) %>% 
  left_join(locations) %>% 
  select(location_name, sum_samples, latitude, longitude) %>% 
  distinct(location_name, .keep_all = TRUE) %>% 
  drop_na()

spatDat <- st_as_sf(spatDat, coords = c("longitude", "latitude"), crs = 4326)

spatDat %>% 
  filter(sum_samples > 15) %>% 
  ggplot() +
  geom_sf(color = "red", size = 2) +
  borders("state", fill = NA, colour = "gray50") +  # US state outlines
  coord_sf(
    xlim = c(-90, -82),
    ylim = c(41.6, 48.3),
    expand = FALSE) +
  theme_minimal() +
  labs(title = "Cisco samples @ MSU n > 15", x = "longitude", y = "latitude")

# View the icollections# View the imported data
head(data_from_access_dbi)

# Close the connection
dbDisconnect(con)