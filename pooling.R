suppressPackageStartupMessages({
  library(tidyverse)
  library(ggplot2)
  library(scales)
  library(patchwork)
  library(vctrs)
  library(xml2)
  library(ggplot2)
  library(terra)
  library(sf)
  library(tmap)
  library(tidycensus)
  library(data.table)
  library(ggspatial)
  library(future)
  library(future.apply)
  library(dplyr)
  library(purrr)
  library(tidyr)
  library(fst)
  library(progressr)
  library(reasdr)
})


datadir <- Sys.getenv("DATA_DIR", "L:/Wetland Flood Mitigation")
pth <- function(...) file.path(datadir, ...)

#-----------------------------------------------------#
# Load banks and service areas
#-----------------------------------------------------#
setwd("L:/Wetland Flood Mitigation")

sas <- readRDS("Service Areas\\ServiceAreas_agg.rds")
banks=readRDS("Bank Footprints/footprints_and_buffers.rds")

banks <- st_as_sf(banks)
sas <- st_as_sf(sas)

banks = banks %>%
  mutate(Name = recode(Name, 
                       "FP_AndL_Everglades_Phase_II_MB" = "FP_L_Everglades_Phase_II_MB"))

sanitize_name <- function(name) {
  name <- gsub("[^A-Za-z0-9_-]", "_", name)
  name <- gsub("_+", "_", name)
  name <- gsub("^_|_$", "", name)
  return(name)
}

summary_bank = readRDS("Extractions and Summaries/Bank Extract/bank_summary.rds")

banks = banks %>% 
  select(c("Name", "geometry"))
sas <- sas %>% 
  select(c("ID", "geometry"))

sas <- sas %>%
  filter(ID %in% summary_bank$bank_name)
banks <- banks %>%
  filter(Name %in% summary_bank$bank_name)


# FIND IDENTICAL SAS
library(purrr)
sas_bbox <- sas %>%
  st_make_valid() %>%
  mutate(
    bbox = map(geometry, ~st_bbox(.)),
    xmin = map_dbl(bbox, ~.[1]),
    ymin = map_dbl(bbox, ~.[2]),
    xmax = map_dbl(bbox, ~.[3]),
    ymax = map_dbl(bbox, ~.[4]),
    # Round to nearest 10 units instead of 0.1
    xmin_r = round(xmin, -1),  # -1 rounds to nearest 10
    ymin_r = round(ymin, -1),
    xmax_r = round(xmax, -1),
    ymax_r = round(ymax, -1)
  )

bbox_duplicates <- sas_bbox %>%
  st_drop_geometry() %>%
  group_by(xmin_r, ymin_r, xmax_r, ymax_r) %>%
  filter(n() > 1) %>%
  summarise(
    identical_sa_ids = paste(ID, collapse = ", "),
    count = n(),
    .groups = "drop"
  )

print(bbox_duplicates)

#######


#-----------------------------------------------------#
# Pooled Bank ECDF Functions for Identical Service Areas
#-----------------------------------------------------#

########

# ------------------------------------------------------------
#               Libraries & helper utilities
# ------------------------------------------------------------
sanitize_name <- function(x) {
  x %>% 
    str_replace_all("[^A-Za-z0-9_-]", "_") %>%   # keep safe chars
    str_replace_all("_+", "_")  %>%              # collapse repeats
    str_remove_all("^_|_$")                      # trim ends
}

compute_ecdf_set <- function(df, vars = c("housing_units",
                                          "housing_value",
                                          "population")) {
  out <- list()
  for (v in vars) {
    if (!v %in% names(df)) next
    vals <- na.omit(df[[v]])
    out[[v]] <- if (length(vals)) ecdf(vals) else NULL
  }
  out
}

# ------------------------------------------------------------
#            1.  pool ⇄ bank cross-walk  (once)
# ------------------------------------------------------------
#  bbox_duplicates is already in memory; if not, readRDS it here
bank_pool_map <- bbox_duplicates %>% 
  mutate(pool_id = paste0("Pool_", row_number())) %>% 
  select(pool_id, identical_sa_ids) %>% 
  separate_rows(identical_sa_ids, sep = ",\\s*") %>% 
  rename(bank_name_raw = identical_sa_ids) %>% 
  mutate(bank_name = sanitize_name(bank_name_raw))

saveRDS(bank_pool_map, "Variables/bank_pool_map.rds")


