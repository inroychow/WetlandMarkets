suppressPackageStartupMessages({
  library(terra)
  library(sf)
  library(dplyr)
  library(future.apply)
  library(progressr)
  library(data.table)
  library(tidyverse)
})

# ------------------------------------------------------------------
#                       Base path & helper
# ------------------------------------------------------------------
datadir <- "C:/Users/indumati/Box/Paper2_final"   
pth     <- function(...) file.path(datadir, ...)  

# ------------------------------------------------------------------
#            Example raster & bank footprints to re‑project
# ------------------------------------------------------------------
example_rast <- rast(pth("Service Area Mosaics", "MT-Running_Colter_Ranch_mosaic.tif"))

banks <- vect(pth("Bank Footprints", "footprints_and_buffers.rds"))
banks <- project(banks, crs(example_rast))   # use same CRS as example raster

# ------------------------------------------------------------------
#                Summarise *bank* extractions
# ------------------------------------------------------------------
bank_extract <- readRDS(pth("Extractions and Summaries",
                            "Extract Banks",
                            "all_bank_extractions.rds"))

bank_summary <- bank_extract %>%
  group_by(bank_name) %>%
  summarise(
    sum_housing_units    = sum(`housing units`, na.rm = TRUE),
    mean_housing_units   = mean(`housing units`, na.rm = TRUE),
    median_housing_units = median(`housing units`, na.rm = TRUE),
    min_housing_units    = min(`housing units`, na.rm = TRUE),
    max_housing_units    = max(`housing units`, na.rm = TRUE),
    sd_housing_units     = sd (`housing units`, na.rm = TRUE),
    q25_housing_units    = quantile(`housing units`, .25, na.rm = TRUE),
    q75_housing_units    = quantile(`housing units`, .75, na.rm = TRUE),
    
    sum_housing_value    = sum(`housing value`, na.rm = TRUE),
    mean_housing_value   = mean(`housing value`, na.rm = TRUE),
    median_housing_value = median(`housing value`, na.rm = TRUE),
    min_housing_value    = min(`housing value`, na.rm = TRUE),
    max_housing_value    = max(`housing value`, na.rm = TRUE),
    sd_housing_value     = sd (`housing value`, na.rm = TRUE),
    q25_housing_value    = quantile(`housing value`, .25, na.rm = TRUE),
    q75_housing_value    = quantile(`housing value`, .75, na.rm = TRUE),
    
    sum_population       = sum(population, na.rm = TRUE),
    mean_population      = mean(population, na.rm = TRUE),
    median_population    = median(population, na.rm = TRUE),
    min_population       = min(population, na.rm = TRUE),
    max_population       = max(population, na.rm = TRUE),
    sd_population        = sd (population, na.rm = TRUE),
    q25_population       = quantile(population, .25, na.rm = TRUE),
    q75_population       = quantile(population, .75, na.rm = TRUE),
    
    cell_count           = n(),
    .groups = "drop"
  )

saveRDS(bank_summary,
        pth("Extractions and Summaries", "bank_summary.rds"))

# ------------------------------------------------------------------
#                Summarise *service‑area* extractions
# ------------------------------------------------------------------
extracted_dir <- pth("Extractions and Summaries", "Extract Service Areas")
output_dir_sa <- pth("Extractions and Summaries", "Summaries Service Areas")
dir.create(output_dir_sa, recursive = TRUE, showWarnings = FALSE)

sa_files <- list.files(extracted_dir, pattern = "^sa_extract_.*\\.rds$", full.names = TRUE)
done     <- list.files(output_dir_sa, pattern = "^summary_.*\\.rds$") |>
  sub("^summary_(.*)\\.rds$", "\\1", x = _)

summarise_sa <- function(df, sa_name) {
  bank_cells <- bank_extract$cell[bank_extract$bank_name == sa_name]
  df <- df[!df$cell %in% bank_cells, ]
  
  df %>%
    summarise(
      sum_housing_units    = sum(housing_units, na.rm = TRUE),
      mean_housing_units   = mean(housing_units, na.rm = TRUE),
      median_housing_units = median(housing_units, na.rm = TRUE),
      min_housing_units    = min(housing_units, na.rm = TRUE),
      max_housing_units    = max(housing_units, na.rm = TRUE),
      sd_housing_units     = sd (housing_units, na.rm = TRUE),
      q25_housing_units    = quantile(housing_units, .25, na.rm = TRUE),
      q75_housing_units    = quantile(housing_units, .75, na.rm = TRUE),
      
      sum_housing_value    = sum(housing_value, na.rm = TRUE),
      mean_housing_value   = mean(housing_value, na.rm = TRUE),
      median_housing_value = median(housing_value, na.rm = TRUE),
      min_housing_value    = min(housing_value, na.rm = TRUE),
      max_housing_value    = max(housing_value, na.rm = TRUE),
      sd_housing_value     = sd (housing_value, na.rm = TRUE),
      q25_housing_value    = quantile(housing_value, .25, na.rm = TRUE),
      q75_housing_value    = quantile(housing_value, .75, na.rm = TRUE),
      
      sum_population       = sum(population, na.rm = TRUE),
      mean_population      = mean(population, na.rm = TRUE),
      median_population    = median(population, na.rm = TRUE),
      min_population       = min(population, na.rm = TRUE),
      max_population       = max(population, na.rm = TRUE),
      sd_population        = sd (population, na.rm = TRUE),
      q25_population       = quantile(population, .25, na.rm = TRUE),
      q75_population       = quantile(population, .75, na.rm = TRUE),
      
      cell_count           = n(),
      .groups = "drop"
    ) %>%
    mutate(sa_id = sa_name)
}

for (fp in sa_files) {
  sa_id <- sub("^sa_extract_(.*)\\.rds$", "\\1", basename(fp))
  if (sa_id %in% done) next
  
  df <- readRDS(fp)
  names(df)[2:4] <- c("housing_units", "housing_value", "population")
  sa_sum <- summarise_sa(df, sa_id)
  
  saveRDS(sa_sum, pth(output_dir_sa, paste0("summary_", sa_id, ".rds")))
}

# Combine every SA summary into one table
sa_summary <- list.files(output_dir_sa, pattern = "^summary_.*\\.rds$", full.names = TRUE) |>
  lapply(readRDS) |>
  bind_rows()

saveRDS(sa_summary, pth("Extractions and Summaries", "sa_summary.rds"))

# ------------------------------------------------------------------
#                Summarise Loss extractions
# ------------------------------------------------------------------
extracted_dir <- pth("Extractions and Summaries", "Extract Loss")
output_dir_loss <- pth("Extractions and Summaries", "Summaries Loss")
dir.create(output_dir_sa, recursive = TRUE, showWarnings = FALSE)

sa_files <- list.files(extracted_dir, pattern = "^loss_extract_.*\\.rds$", full.names = TRUE)
done     <- list.files(output_dir_loss, pattern = "^summaryloss_.*\\.rds$") |>
  sub("^summary_(.*)\\.rds$", "\\1", x = _)

summarise_loss <- function(df, sa_name) {
  bank_cells <- bank_extract$cell[bank_extract$bank_name == sa_name]
  df <- df[!df$cell %in% bank_cells, ]
  
  df %>%
    summarise(
      sum_housing_units    = sum(housing_units, na.rm = TRUE),
      mean_housing_units   = mean(housing_units, na.rm = TRUE),
      median_housing_units = median(housing_units, na.rm = TRUE),
      min_housing_units    = min(housing_units, na.rm = TRUE),
      max_housing_units    = max(housing_units, na.rm = TRUE),
      sd_housing_units     = sd (housing_units, na.rm = TRUE),
      q25_housing_units    = quantile(housing_units, .25, na.rm = TRUE),
      q75_housing_units    = quantile(housing_units, .75, na.rm = TRUE),
      
      sum_housing_value    = sum(housing_value, na.rm = TRUE),
      mean_housing_value   = mean(housing_value, na.rm = TRUE),
      median_housing_value = median(housing_value, na.rm = TRUE),
      min_housing_value    = min(housing_value, na.rm = TRUE),
      max_housing_value    = max(housing_value, na.rm = TRUE),
      sd_housing_value     = sd (housing_value, na.rm = TRUE),
      q25_housing_value    = quantile(housing_value, .25, na.rm = TRUE),
      q75_housing_value    = quantile(housing_value, .75, na.rm = TRUE),
      
      sum_population       = sum(population, na.rm = TRUE),
      mean_population      = mean(population, na.rm = TRUE),
      median_population    = median(population, na.rm = TRUE),
      min_population       = min(population, na.rm = TRUE),
      max_population       = max(population, na.rm = TRUE),
      sd_population        = sd (population, na.rm = TRUE),
      q25_population       = quantile(population, .25, na.rm = TRUE),
      q75_population       = quantile(population, .75, na.rm = TRUE),
      
      cell_count           = n(),
      .groups = "drop"
    ) %>%
    mutate(sa_id = sa_name)
}

for (fp in sa_files) {
  sa_id <- sub("^loss_extract_(.*)\\.rds$", "\\1", basename(fp))
  if (sa_id %in% done) next
  
  df <- readRDS(fp)
  names(df)[2:4] <- c("housing_units", "housing_value", "population")
  sa_sum <- summarise_sa(df, sa_id)
  
  saveRDS(sa_sum, pth(output_dir_sa, paste0("summaryloss_", sa_id, ".rds")))
}

# Combine every SA summary into one table
loss_summary <- list.files(output_dir_sa, pattern = "^summaryloss_.*\\.rds$", full.names = TRUE) |>
  lapply(readRDS) |>
  bind_rows()

saveRDS(loss_summary, pth("Extractions and Summaries", "loss_summary.rds"))

# ------------------------------------------------------------------
#                  Combine bank, SA, loss summaries
# ------------------------------------------------------------------

bank_tbl <- bank_summary  |>
  rename(id = bank_name,
         sum_hu_bank    = sum_housing_units,
         mean_hu_bank   = mean_housing_units,
         median_hu_bank = median_housing_units,
         max_hu_bank    = max_housing_units) |>
  select(id, starts_with("sum_"), starts_with("mean_"),
         starts_with("median_"), starts_with("max_"))

sa_tbl <- sa_summary |>
  rename(id = sa_id,
         sum_hu_sa    = sum_housing_units,
         mean_hu_sa   = mean_housing_units,
         median_hu_sa = median_housing_units,
         max_hu_sa    = max_housing_units) |>
  select(id, starts_with("sum_"), starts_with("mean_"),
         starts_with("median_"), starts_with("max_"))


full_summary <- bank_tbl %>%
  inner_join(sa_tbl,   by = "id") %>%

saveRDS(full_summary,
        pth("Extractions and Summaries", "full_summary.rds"))
