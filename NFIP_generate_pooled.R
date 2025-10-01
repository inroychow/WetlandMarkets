setwd("C:/Users/indumati/Box/Paper2_final")

library(pacman)
p_load(tidyverse, rgdal, sp, raster, sf, RColorBrewer, ggplot2, dplyr, raster, gstat, dismo, openxlsx, lubridate, data.table)

# ------------------------------- Load all necessary raw data -----------------------------

# NFIP data: FIMA NFIP Redacted Claims/Policies Dataset
nfipclaims <- read_csv("NFIP/Claims/FimaNfipClaims.csv")
nfippolicies <- fread("NFIP/Policies/FimaNfipPolicies.csv")

# Bank and huc12 geometry
banks = readRDS("Footprints/footprints_and_buffers.rds")
huc12s = readRDS("CONUS/Variables/huc12_boundaries.rds")
huc12s = st_transform(huc12s, st_crs(banks))

huc12s = st_as_sf(huc12s)
banks = st_as_sf(banks)

banks <- st_make_valid(banks)
huc12s <- st_make_valid(huc12s)

# HUC12 to downstream tracts mapping 
htdt = readRDS("CONUS/Variables/htdt_conus.rds")

# ------------------------------- CPI to convert to 2025 dollars ----------------------------------

cpi_df <- tibble(
  year = 1978:2025,
  cpi = c(65.2, 72.6, 82.4, 90.9, 96.5, 99.6, 103.9, 107.6, 109.6, 113.6, 118.3,
          124.0, 130.7, 136.2, 140.3, 144.5, 148.2, 152.4, 156.9, 160.5, 163.0,
          166.6, 172.2, 177.1, 179.9, 184.0, 188.9, 195.3, 201.6, 207.3, 215.3,
          214.5, 218.1, 224.9, 229.6, 232.957, 236.736, 237.017, 240.007, 245.120,
          251.107, 255.657, 258.811, 270.970, 292.655, 301.808, 310.350, 318.000)
)

# ---------------------------------------------------
# ----------------- CLAIMS --------------------------
# ---------------------------------------------------

claims <- nfipclaims[,c("censusTract", "yearOfLoss",
                        "amountPaidOnBuildingClaim", "amountPaidOnContentsClaim", 
                        "totalBuildingInsuranceCoverage", "totalContentsInsuranceCoverage")]

colnames(claims) <- c("census_tract", "year", "buildings_paid", "contents_paid", "buildings_coverage", "contents_coverage")
claims$claims_count <- 1

claims = claims %>% 
  mutate(total_paid = buildings_paid + contents_paid)

# Adjust for inflation
claims_adj <- claims %>%
  left_join(cpi_df, by = "year") %>%
  mutate(cpi_2025 = 318.000,
         inflation_factor = cpi_2025 / cpi,
         total_paid_2025 = total_paid * inflation_factor)

# Round to nearest dollar
claims_adj$total_paid_2025 <- round(claims_adj$total_paid_2025)

# Unnest downstream tracts from htdt
tracts_to_keep <- htdt %>%
  unnest(cols = downstream_tracts) %>%
  distinct(downstream_tracts) %>%
  pull(downstream_tracts)

# Subset claims_adj to only include those census tracts
claims <- claims_adj %>%
  filter(census_tract %in% tracts_to_keep) %>% 
  rename(tract = census_tract)

claims = claims %>% 
  dplyr::select(-cpi, -cpi_2025, -inflation_factor)

# Perform intersection: returns geometry of overlaps
intersection <- st_intersection(
  banks %>% dplyr::select(Name),
  huc12s %>% dplyr::select(huc12)
)

bank_huc12 <- intersection %>%
  st_drop_geometry() %>%
  distinct(Name, huc12)

# Expand the downstream tract list into one row per pair
htdt_long <- htdt %>%
  tidyr::unnest_longer(downstream_tracts) %>%
  rename(tract = downstream_tracts)

# ------------------------------------------------------------------
# Set up projected CRS for distance calculations
# ------------------------------------------------------------------
proj_crs <- 5070  # US National Albers (metres)

# ------------------------------------------------------------------
# Prepare geometries
# ------------------------------------------------------------------
banks_proj <- st_transform(banks, proj_crs)
tracts_centroids <- st_read("CONUS/Variables/tracts_centroids.shp", quiet = TRUE) %>%
  st_transform(proj_crs)

# ------------------------------------------------------------------
# Build the initial tract-bank table (hydrologically upstream)
# ------------------------------------------------------------------
banks_upstream_tracts <- bank_huc12 %>%
  inner_join(htdt_long, by = "huc12") %>%
  distinct(tract, Name)

# ------------------------------------------------------------------
# Add geometries and compute straight-line distance (km)
# ------------------------------------------------------------------
banks_upstream_geom <- banks_upstream_tracts %>%
  left_join(banks_proj %>% select(Name, bank_geom = geometry),
            by = "Name") %>%
  left_join(tracts_centroids %>% select(tract = GEOID, tract_geom = geometry),
            by = "tract") %>%
  mutate(dist_km = as.numeric(
    st_distance(bank_geom, tract_geom, by_element = TRUE)) / 1000)

# ------------------------------------------------------------------
# Keep only links where the bank lies ≤ 40 km upstream
# ------------------------------------------------------------------
banks_upstream_tracts_with_dist <- banks_upstream_geom %>%
  filter(dist_km <= 40) %>%
  st_drop_geometry() %>%
  select(tract, Name, dist_km)  # KEEP dist_km here
# ------------------------------------------------------------------
# Create claims by tract-year
# ------------------------------------------------------------------
claims_by_tract_year <- claims %>%
  group_by(tract, year) %>%
  summarise(
    n_claims = n(), 
    total_paid_2025 = sum(total_paid_2025, na.rm = TRUE),
    .groups = "drop"
  )

# ---------------------------------------------------
# ----------------- POLICIES ------------------------
# ---------------------------------------------------
policies_x = nfippolicies %>% 
  rename(policies_count = policyCount,
         buildings_coverage = totalBuildingInsuranceCoverage,
         contents_coverage = totalContentsInsuranceCoverage,
         tract = censusTract) %>% 
  mutate(total_coverage = buildings_coverage + contents_coverage) %>% 
  mutate(policyeffective_year = year(ymd(policyEffectiveDate))) %>% 
  dplyr::select(tract, policies_count, policyeffective_year, buildings_coverage, contents_coverage, total_coverage)

# Join to claims data and adjust for inflation
policies_adj <- policies_x %>%
  left_join(cpi_df, by = c("policyeffective_year" = "year")) %>%
  mutate(cpi_2025 = 318.000,
         inflation_factor = cpi_2025 / cpi,
         total_coverage_2025 = total_coverage * inflation_factor)

policies_adj$total_coverage_2025 <- round(policies_adj$total_coverage_2025)

policies_adj$tract = str_pad(policies_adj$tract, pad = "0", width = 11)
policies_summary <- policies_adj %>%
  group_by(tract, policyeffective_year) %>%
  summarize(
    policies_count = sum(policies_count, na.rm = TRUE),
    buildings_coverage = sum(buildings_coverage, na.rm = TRUE),
    contents_coverage = sum(contents_coverage, na.rm = TRUE),
    total_coverage_2025 = sum(total_coverage_2025, na.rm = TRUE),
    .groups = "drop"
  )

# Save summarized data
saveRDS(policies_summary, "NFIP/Policies/Policies_summary.rds")

# ------------------------------------------------------------------
# POOLED VERSION WITH DISTANCE METRICS
# ------------------------------------------------------------------

# Read the pool map
poolmap <- readRDS("CONUS/HUC8 Validation/pool_map.rds")

# Add pool_id to tract-bank pairs and preserve distance information
banks_upstream_tracts_pooled <- banks_upstream_tracts_with_dist %>%
  left_join(poolmap %>% select(pool_id, bank_name),
            by = c("Name" = "bank_name")) %>%
  mutate(pool_id = coalesce(pool_id, paste0("Solo_", Name)))
# First, get total pool sizes
pool_sizes <- poolmap %>%
  group_by(pool_id) %>%
  summarise(total_banks_in_pool = n_distinct(bank_name), .groups = "drop")

# Then join to your pools_per_tract
pools_per_tract <- banks_upstream_tracts_pooled %>%
  group_by(tract, pool_id) %>%
  summarise(
    banks_in_pool = paste(unique(Name), collapse = "; "),
    n_banks_upstream = n_distinct(Name),  # renamed for clarity
    min_dist_km = min(dist_km),
    mean_dist_km = mean(dist_km),
    max_dist_km = max(dist_km),
    .groups = "drop"
  ) %>%
  left_join(pool_sizes, by = "pool_id") %>%
  mutate(total_banks_in_pool = coalesce(total_banks_in_pool, 1))  # solo banks = 1


# Create final claims dataset with pools
claims_final_pooled <- claims_by_tract_year %>%
  left_join(pools_per_tract, by = "tract") %>%
  left_join(policies_summary, by = c("tract", "year" = "policyeffective_year")) %>%
  rename(claim_year = year,
         upstream_bank_or_pool = pool_id) %>%
  mutate(upstream_bank_or_pool = str_remove(upstream_bank_or_pool, "^Solo_")) %>%
  filter(!is.na(upstream_bank_or_pool)) %>%
  select(-buildings_coverage, -contents_coverage)

claims_final = claims_final_pooled %>% 
  select(tract, claim_year, n_claims, total_paid_2025, upstream_bank_or_pool, n_banks_upstream, policies_count, total_coverage_2025)

# Save final dataset
saveRDS(claims_final, "NFIP/nfipbankdata_pool.rds")

