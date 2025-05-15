#We are subsetting and generating the NFIP claims data on a county level and a zipcode level.

setwd("C:/Users/indumati/Box/Mitigation bank generate backup 5.9")
library(pacman)
p_load(tidyverse, rgdal, sp, raster, sf, RColorBrewer, ggplot2, dplyr, raster, gstat, dismo, openxlsx)

#Load NFIP claims CSV - Source: FIMA NFIP Redacted Claims Dataset
nfip <- read_csv("NFIP/Claims/FimaNfipClaims.csv")

claims <- nfip[,c("countyCode", "reportedZipCode", "yearOfLoss", "dateOfLoss",
                    "amountPaidOnBuildingClaim", "amountPaidOnContentsClaim", 
                    "totalBuildingInsuranceCoverage", "totalContentsInsuranceCoverage", "crsClassificationCode", "ratedFloodZone", "waterDepth", "latitude", "longitude", "censusTract")]

colnames(claims) <- c("county_fip", "zip", "year", "date", "buildings_paid", "contents_paid", "buildings_coverage", "contents_coverage", "CRS", "floodZone", "water_depth", "lat", "long", "census_tract")
claims$claims_count <- 1

#Deal with CRS
claims$CRS <- as.character(claims$CRS)
claims$floodZone <- trimws(claims$floodZone)
# Re-classify flood zones into SFHA and Non-SFHA
SFHA <- c("AE", "VE", "V", "A", "AOB", "A14", "AO", "A07", "A09", "A05", "A08", 
          "A06", "A04", "A10", "A13", "A01", "AH", "A11", "A17", "AHB", "A03", "V20", "A15", "A12", "V15", 
          "D", "A02", "V08", "A16", "V09", "V10", "V19", "V12", "V11", "A30", "A18", "A21", "A22", "V13", 
          "V16", "A20", "A19", "V06", "V18", "V17", "A24", "A25", "V05", "V14", "V04", "V21", "V07", "V02", "A23", 
          "A26", "A28", "V03", "V01", "V22", "A27", "V23", "V30", "V24", "V27")

non_SFHA <- c("B", "C", "X", "AR", "A99") #AR And A99 are not treated as SFHA for CRS classification
claims$floodZone <- trimws(claims$floodZone)
claims$floodZone <- ifelse(claims$floodZone %in% SFHA, "SFHA", 
                      ifelse(claims$floodZone %in% non_SFHA, "Non SFHA", "Unknown"))

crs_discount_SFHA <- c('1' = 0.45, '2' = 0.40, '3' = 0.35, '4' = 0.30, '5' = 0.25, '6' = 0.20, '7' = 0.15, '8' = 0.10, '9' = 0.05, '10' = 0.00)

crs_discount_nonSFHA <- c('1' = 0.10, '2' = 0.10, '3' = 0.10, '4' = 0.10, '5' = 0.10, '6' = 0.10, '7' = 0.05, '8' = 0.05, '9' = 0.05, '10' = 0.00)

claims$CRS_discount <- ifelse(claims$floodZone == "SFHA", 
                              crs_discount_SFHA[claims$CRS],
                              ifelse(claims$floodZone == "Non SFHA",
                                     crs_discount_nonSFHA[claims$CRS], NA))

claims$CRS_discount[is.na(claims$CRS_discount)] <- 0

claims = claims %>% 
  mutate(total_paid = buildings_paid + contents_paid)
saveRDS(claims, "NFIP/Claims/claims.rds")




# Do florida analysis

library(tigris)  # to get census tract shapefiles

# Load Florida tract shapefile
fl_tracts <- tracts(state = "FL", year = 2020, class = "sf")
fl_tracts = st_transform(fl_tracts, 4326)
# Convert your claims data to an sf object using coordinates
claims_sf <- claims %>%
  filter(!is.na(lat) & !is.na(long)) %>%
  st_as_sf(coords = c("long", "lat"), crs = 4326, remove = FALSE)

# Spatial join to get census tracts from FL only
claims_fl <- st_join(claims_sf, fl_tracts["GEOID"])



# Filter for 1983–1987 and 2019–2023
claims_window <- claims_fl %>%
  filter(year %in% c(1983:1987, 2019:2023))

# Create a column for the time window (1985 or 2021)
claims_window <- claims_window %>%
  mutate(window = case_when(
    year %in% 1983:1987 ~ "1985_5yr",
    year %in% 2019:2023 ~ "2021_5yr",
    TRUE ~ NA_character_
  ))

# Group by census tract and window, then calculate average claim amount
avg <- claims_window %>%
  group_by(GEOID, window) %>%
  summarise(
    avg_total_paid = mean(total_paid, na.rm = TRUE),
    avg_claims_count = mean(claims_count, na.rm = TRUE),
    .groups = "drop"
  )


avg_adj <- avg %>%
  mutate(
    avg_total_paid_2021 = case_when(
      window == "1985_5yr" ~ avg_total_paid * 2.52,
      window == "2021_5yr" ~ avg_total_paid  # already in 2021 dollars
    ),
    avg_claims_count_2021 = avg_claims_count  
  )


library(tigris)
library(dplyr)
library(tidyr)
library(sf)
library(ggplot2)
library(scales)

# --- 0. (Re-)load tract shapes, keep as sf ---------------------------
fl_tracts <- tracts(state = "FL", year = 2020, class = "sf") |>
  st_transform(4326)

# --- 1. Make the wide table & drop geometry --------------------------
diff_df <- avg_adj |>                          # from your previous step
  st_drop_geometry() |>                        # <-- key line
  dplyr::select(GEOID, window, avg_total_paid_2021) |>
  pivot_wider(names_from  = window,
              values_from = avg_total_paid_2021) |>
  mutate(diff_paid_2021 = `2021_5yr` - `1985_5yr`)

# --- 2. Join attributes onto tracts ----------------------------------
fl_tracts_diff <- fl_tracts |>
  left_join(diff_df, by = "GEOID")             # diff_df is now plain data

# --- 3. Plot ----------------------------------------------------------
ggplot(fl_tracts_diff) +
  geom_sf(aes(fill = diff_paid_2021),
          color = NA, linewidth = 0) +
  scale_fill_gradient2(
    name = "Δ Avg. Paid\n(2021 $)",
    low  = muted("blue"), mid = "white", high = muted("red"),
    midpoint = 0,
    labels = label_dollar()
  ) +
  labs(
    title    = "Change in Average NFIP Total Paid per Census Tract",
    subtitle = "Mean of 2019-23 window minus mean of 1983-87 window (2021 dollars)",
    caption  = "Source: NFIP claims; tract shapes: TIGER/Line 2020"
  ) +
  theme_minimal() +
  theme(panel.grid.major = element_blank(),
        axis.title       = element_blank())


### b y zip 
library(tidycensus)
# 1. Florida ZCTA (ZIP) polygons -------------------------------------
options(tigris_use_cache = TRUE)

# 1. Get all ZCTA geometries (nationwide)
zips_all <- zctas(year = 2021)

# 2. Get Florida state boundary
fl_state <- states(cb = TRUE, year = 2020) %>%
  filter(STUSPS == "FL") %>%
  st_transform(st_crs(zips_all))

# 3. Spatial filter: keep only ZCTAs that intersect Florida
fl_zips <- st_filter(zips_all, fl_state)

zip_id = "ZCTA5CE20"

fl_zips = st_transform(fl_zips, st_crs(claims_fl))
# 2. Spatial join: each claim gets its ZIP ---------------------------
claims_fl_zip <- st_join(claims_sf, fl_zips[zip_id]) |>
  filter(!is.na(.data[[zip_id]]))               # drop claims that miss a ZIP

# 3. Define the time windows -----------------------------------------
claims_zip_win <- claims_fl_zip |>
  filter(year %in% c(1983:1987, 2019:2023)) |>
  mutate(window = case_when(
    year %in% 1983:1987 ~ "1985_5yr",
    year %in% 2019:2023 ~ "2021_5yr"
  ))

# 4. Average $ and count per ZIP & window ----------------------------
avg_zip <- claims_zip_win |>
  group_by(.data[[zip_id]], window) |>
  summarise(
    avg_total_paid = mean(total_paid,  na.rm = TRUE),
    avg_claims_cnt = mean(claims_count, na.rm = TRUE),
    .groups = "drop"
  )

# 5. Inflate the older $$ to 2021 dollars ----------------------------
avg_zip_adj <- avg_zip |>
  mutate(avg_total_paid_2021 = if_else(window == "1985_5yr",
                                       avg_total_paid * 2.52,   # CPI factor
                                       avg_total_paid),
         avg_claims_cnt_2021 = avg_claims_cnt)

# 6. Pivot wider & compute difference -------------------------------
diff_zip <- avg_zip_adj |>
  dplyr::select(all_of(zip_id), window, avg_total_paid_2021) |>
  pivot_wider(names_from  = window,
              values_from = avg_total_paid_2021) |>
  mutate(diff_paid_2021 = `2021_5yr` - `1985_5yr`) |>
  st_drop_geometry()                             # ordinary data frame

# 7. Join attributes onto ZCTA shapes --------------------------------
fl_zips_diff <- fl_zips |>
  left_join(diff_zip, by = setNames(zip_id, zip_id))  # keyed on ZIP

# 8. Map --------------------------------------------------------------
ggplot(fl_zips_diff) +
  geom_sf(aes(fill = diff_paid_2021), color = NA, linewidth = 0) +
  scale_fill_gradient2(
    name = "Δ Avg. Paid\n(2021 $)",
    low = muted("blue"), mid = "white", high = muted("red"),
    midpoint = 0, labels = label_dollar()
  ) +
  labs(
    title    = "Change in Average NFIP Total Paid per ZIP (ZCTA)",
    subtitle = "Mean of 2019-23 window minus mean of 1983-87 window (2021 dollars)",
    caption  = "Source: NFIP claims; ZCTA shapes: TIGER/Line 2020"
  ) +
  theme_minimal() +
  theme(panel.grid.major = element_blank(),
        axis.title       = element_blank())
