
setwd("C:/Users/indumati/Box/Paper2_final")


library(pacman)
p_load(tidyverse, rgdal, sp, raster, sf, RColorBrewer, ggplot2, dplyr, raster, gstat, dismo, openxlsx, lubridate, data.table)


# ------------------------------- Load all necessary raw data -----------------------------

# NFIP data: FIMA NFIP Redacted Claims/Policies Dataset -----
nfipclaims <- read_csv("NFIP/Claims/FimaNfipClaims.csv")
nfippolicies <- fread("CONUS/FimaNfipPolicies.csv") # use data.table package bc large dataset

# Bank and huc12 geometry -----
banks = readRDS("Bank Footprints/footprints_and_buffers.rds")
huc12s = readRDS("Variables/huc12_boundaries.rds")
huc12s = st_transform(huc12s, st_crs(banks))

huc12s = st_as_sf(huc12s)
banks = st_as_sf(banks)

banks <- st_make_valid(banks)
huc12s <- st_make_valid(huc12s)

# HUC12 to downstream tracts mapping 
htdt = readRDS("Variables/htdt_conus.rds")

#------------------------------- CPI to convert to 2025 dollars ----------------------------------

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

# -----------------------------------------

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
  dplyr::select(-cpi, -cpi_2025, -lat, -inflation_factor)


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

banks_upstream_tracts <- bank_huc12 %>%
  inner_join(htdt_long, by = "huc12") %>%
  distinct(tract, Name)

claims_by_tract_year <- claims %>%              # start from the claim-level data
  group_by(tract, year) %>%                         # <- one row per tract *and* year
  summarise(
    n_claims        = n(),                          # how many claims
    total_paid_2025 = sum(total_paid_2025,  na.rm = TRUE),
    .groups = "drop"
  )

## ----------------  banks: 1 row per tract ----------------
banks_per_tract <- banks_upstream_tracts %>%        # already tract–bank pairs
  group_by(tract) %>%                               # collapse to one row per tract
  summarise(
    banks_upstream  = paste(unique(Name), collapse = "; "),   # list of names
    n_banks         = n_distinct(Name),
    .groups = "drop"
  )


## ---------------- 3. final join ----------------
claims_final <- claims_by_tract_year %>%
  left_join(banks_per_tract, by = "tract")

sanitize_path <- function(path) {
  path <- gsub("[[:space:]]+", "_", path)
  path <- gsub("[^[:alnum:]_/-]", "", path)
  return(path)
}

# Apply it to semicolon-separated names
claims_final <- claims_final %>%
  mutate(
    banks_upstream = banks_upstream %>%
      str_split(";\\s*") %>%                            # split on semicolons
      map_chr(~ paste(map_chr(.x, sanitize_path),      # sanitize each name
                      collapse = ";"))                 # rejoin with semicolon
  ) %>% 
  rename(claim_year = year)


# ---------------------------------------------------
# ----------------- POLICIES ------------------------
# ---------------------------------------------------
policies_x = policies %>% 
  rename(policies_count = policyCount,
         buildings_coverage = totalBuildingInsuranceCoverage,
         contents_coverage = totalContentsInsuranceCoverage,
         tract = censusTract) %>% 
  mutate(total_coverage = buildings_coverage + contents_coverage) %>% 
  mutate(policyeffective_year = year(ymd(policyEffectiveDate))) %>% 
  dplyr::select(tract, policies_count, policyeffective_year, buildings_coverage, contents_coverage, total_coverage)



# Join to claims data and adjust for inflation
policies_adj <- policies_x %>%
  left_join(cpi_df, by = "policyeffective_year") %>%
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
saveRDS(policies_summary, "NFIP\\Policies")


claims_final = claims_final %>% 
  left_join(policies_summary, by = c("tract", "claim_year" = "policyeffective_year"))

claims_final <- claims_final %>%
  tidyr::separate_rows(banks_upstream, sep = ";") 

claims_final = claims_final %>% 
  dplyr::select(-c(buildings_coverage, contents_coverage))

saveRDS(claims_final, "NFIP/nfipbankdata.rds")

