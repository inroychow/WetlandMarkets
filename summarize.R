library(terra)
library(progressr)
library(sf)
library(dplyr)
library(future.apply)

setwd("C:/Users/indumati/Box/Paper2_final")

# ---------------------------------------------- #
# ------------- Summarize banks -----------------# 
# ---------------------------------------------- #

example_rast = rast("Mosaics/MT-Running_Colter_Ranch_mosaic.tif")

banks = vect("Bank Footprints/footprints_and_buffers.rds")
banks = project(banks, crs(final_multilayer))

# --------------------------------------------------

# Paths to spatial datasets
huc12_shp <- "Variables/huc12_boundaries.gpkg"
sa_shp <- "Service Areas/ServiceAreas_agg.rds"
output_folder <- "Service Areas/Temp lookup"


## Summarize bank -----------------------------------------------------------------------
bank_extract = readRDS("bank_cellvalues.rds")
# colnames(bank_extract)[2:4] <- c("housing units", "housing value", "population")
bank_summary <- bank_extract %>%
  # Exclude rows with no raster data
  group_by(bank_name) %>%
  summarise(
    sum_housing_units = sum(`housing units`, na.rm = TRUE),
    mean_housing_units = mean(`housing units`, na.rm = TRUE),
    median_housing_units = median(`housing units`, na.rm = TRUE),
    min_housing_units = min(`housing units`, na.rm = TRUE),
    max_housing_units = max(`housing units`, na.rm = TRUE),
    sd_housing_units = sd(`housing units`, na.rm = TRUE),
    q25_housing_units = quantile(`housing units`, 0.25, na.rm = TRUE),
    q75_housing_units = quantile(`housing units`, 0.75, na.rm = TRUE),
    
    sum_housing_value = sum(`housing value`, na.rm = TRUE),
    mean_housing_value = mean(`housing value`, na.rm = TRUE),
    median_housing_value = median(`housing value`, na.rm = TRUE),
    min_housing_value = min(`housing value`, na.rm = TRUE),
    max_housing_value = max(`housing value`, na.rm = TRUE),
    sd_housing_value = sd(`housing value`, na.rm = TRUE),
    q25_housing_value = quantile(`housing value`, 0.25, na.rm = TRUE),
    q75_housing_value = quantile(`housing value`, 0.75, na.rm = TRUE),
    
    sum_population = sum(population, na.rm = TRUE),
    mean_population = mean(population, na.rm = TRUE),
    median_population = median(population, na.rm = TRUE),
    min_population = min(population, na.rm = TRUE),
    max_population = max(population, na.rm = TRUE),
    sd_population = sd(population, na.rm = TRUE),
    q25_population = quantile(population, 0.25, na.rm = TRUE),
    q75_population = quantile(population, 0.75, na.rm = TRUE),
    
    cell_count = n(),  # Total number of cells
    .groups = "drop"
  )


# saveRDS(bank_summary, "bank_summary.rds")
saveRDS(bank_summary, "bank_summary.rds")



# ---------------------------------------------- #
# ------------- Summarize SAs - -----------------# 
# ---------------------------------------------- #

SAs = readRDS("Service Areas/ServiceAreas_agg.rds")



