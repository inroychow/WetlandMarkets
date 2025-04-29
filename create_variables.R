library(tidyverse)
library(xml2)
library(ggplot2)
library(terra)
library(sf)
library(tmap)
library(tidycensus)
library(data.table)

setwd("C:/Users/indumati/Box/Paper2_final")

# Output folder
main_folder <- "Variables/State-by-state"
dir.create(main_folder, recursive = TRUE, showWarnings = FALSE)

# Load wetland raster
wetland_raster <- rast("LCMAP Rasters/LCMAP_CU_2021_V13_LCPRI.tif")

# Load HUC12 boundaries
huc12_boundaries <- st_read("Geo Boundaries/HUC12s/WBD_National_GDB.gdb", layer = "WBDHU12")
huc12_boundaries <- st_transform(huc12_boundaries, crs(wetland_raster))
huc12_boundaries <- st_make_valid(huc12_boundaries)
# Load CONUS boundaries

conus = readRDS("Geo Boundaries/USA States/conus.rds")
conus <- project(conus, crs(wetland_raster))

# Load service area shapefiles
serviceareas <- vect("Service Areas/ServiceAreasShapefile.shp")
serviceareas <- project(serviceareas, crs(wetland_raster))

# Read upstream data
upstream_df <- fread("Python/service_area_huc12_to_downstream_tracts_2021.csv")

# List of states to process
states = unique(conus$STATEFP)

for (state in states) {
  
  message("Processing state: ", state)
  
  # Get state boundary
  state_boundary <- conus[conus$STATEFP == state, ]  
  
  # Subset service areas to state
  state_SAs <- mask(serviceareas, state_boundary)
  state_SAs$service_area <- gsub("^Document/", "", state_SAs$FolderPath)
  
  # Crop wetland raster to state boundary
  state_wetlands <- crop(wetland_raster, state_boundary)
  
  # Get tract-level data for the state
  hv_tracts <- get_acs(
    geography = "tract", 
    variables = "B25077_001",  # Median home value
    state = state,             
    year = 2021,               
    geometry = TRUE
  ) %>%
    rename(property_value = estimate) %>%
    st_transform(crs(state_wetlands))
  
  hu_tracts <- get_acs(
    geography = "tract", 
    variables = "B25001_001",  # Number of housing units
    state = state,             
    year = 2021,               
    geometry = TRUE
  ) %>%
    rename(housing_units = estimate) %>%
    st_transform(crs(state_wetlands))
  
  pop_tracts <-
    get_acs(
      geography = "tract", 
      variables = "B01003_001",  # Number of housing units
      state = state,             
      year = 2021,               
      geometry = TRUE
    ) %>%
    rename(population = estimate) %>%
    st_transform(crs(state_wetlands))
  
  # Filter HUC12 boundaries to state extent
  state_huc12_boundaries <- st_crop(huc12_boundaries, state_boundary)
  state_huc12_boundaries <- st_transform(state_huc12_boundaries, crs(state_wetlands))
  
  # Rasterize HUC12 boundaries
  state_huc12_vect <- vect(state_huc12_boundaries)
  state_huc12_raster <- rasterize(
    state_huc12_vect,
    state_wetlands,
    field = "huc12",
    touches = TRUE
  )
  
  # Mask HUC12 raster with wetland raster
  state_wetland_huc12_raster <- mask(state_huc12_raster, state_wetlands)
  
  # Calculate tract centroids
  tracts_centroids <- st_centroid(hu_tracts) %>%
    select(GEOID, housing_units) %>%
    st_as_sf()
  
  # Improved logic for state_huc12_to_downstream_tracts
  upstream_df = upstream_df %>% 
    mutate(huc12 = str_pad(huc12, width = 12, pad = "0"))
  
  upstream_subset <- upstream_df %>%
    filter(huc12 %in% state_huc12_boundaries$huc12)
  
  state_huc12_to_downstream_tracts <- upstream_subset %>%
    mutate(
      downstream_tracts = str_remove_all(downstream_tracts, "[\\[\\]']") %>% # Remove brackets and quotes
        str_split(",\\s*")  # Split into lists
    ) %>%
    group_by(huc12) %>%
    summarize(downstream_tracts = list(unique(unlist(downstream_tracts))), .groups = "drop") %>% # Flatten and deduplicate
    mutate(
      huc12 = str_pad(huc12, width = 12, pad = "0"), # Pad HUC12 codes
      downstream_tracts = lapply(downstream_tracts, function(x) str_pad(x, width = 11, pad = "0")) # Pad tract codes
    )
  
  # Save outputs for the state
  state_folder <- paste0(main_folder, "/", state, "/")
  dir.create(state_folder, recursive = TRUE, showWarnings = FALSE)
  
  writeRaster(state_wetlands, paste0(state_folder, "wetlands.tif"), overwrite = TRUE)
  writeRaster(state_huc12_raster, paste0(state_folder, "huc12_raster.tif"), overwrite = TRUE)
  st_write(tracts_centroids, paste0(state_folder, "tracts_centroids.shp"), delete_dsn = TRUE)
  st_write(state_huc12_boundaries, paste0(state_folder, "huc12_boundaries.gpkg"), append = FALSE)
  saveRDS(state_huc12_to_downstream_tracts, paste0(state_folder, "huc12_to_downstream_tracts_2021.rds"))
  saveRDS(hu_tracts, paste0(state_folder, "hu_tracts"))
  saveRDS(hv_tracts, paste0(state_folder, "hv_tracts"))
  saveRDS(pop_tracts, paste0(state_folder, "pop_tracts"))
  
  message("Completed processing for state: ", state)
}

#########
# COMBINING THE VARIABLE DATA FOR STATES, SINCE RELATIONSHIPS CROSS STATE BOUNDARIES 
#########

# Example list of states to process
states = unique(conus$STATEFP)

# Combine all `huc12_to_downstream_tracts` files into one
all_huc12_to_downstream_tracts <- do.call(
  bind_rows,
  lapply(states, function(state) {
    readRDS(paste0("Variables/State-by-state", state, "/huc12_to_downstream_tracts_2021.rds"))
  })
)

# Ensure unique entries
all_huc12_to_downstream_tracts <- distinct(all_huc12_to_downstream_tracts)

# Combine tract centroids for all states
all_tracts_centroids <- do.call(
  bind_rows,
  lapply(states, function(state) {
    st_read(paste0("Variables - all states/", state, "/tracts_centroids.shp"))
  })
)

# Combine attribute data (e.g., housing units, property value, population)
all_hu_tracts <- do.call(
  bind_rows,
  lapply(states, function(state) {
    readRDS(paste0("Variables - all states/", state, "/hu_tracts"))
  })
)

all_hv_tracts <- do.call(
  bind_rows,
  lapply(states, function(state) {
    readRDS(paste0("Variables - all states/", state, "/hv_tracts"))
  })
)

all_pop_tracts <- do.call(
  bind_rows,
  lapply(states, function(state) {
    readRDS(paste0("Variables - all states/", state, "/pop_tracts"))
  })
)


##### -----------------------------------------------HUC12 to downstream tracts CONUS ### -------------------------------------------------


###### Filter out huc12s that are not CONUS ####
# htdt from Python

huc12_boundaries = st_read("Variables/huc12_boundaries.gpkg", quiet = TRUE)

# Define the CONUS state codes
conus_states <- c("AL", "AZ", "AR", "CA", "CO", "CT", "DE", "FL", "GA", "ID", "IL", "IN", "IA",
                  "KS", "KY", "LA", "ME", "MD", "MA", "MI", "MN", "MS", "MO", "MT", "NE", "NV",
                  "NH", "NJ", "NM", "NY", "NC", "ND", "OH", "OK", "OR", "PA", "RI", "SC", "SD",
                  "TN", "TX", "UT", "VT", "VA", "WA", "WV", "WI", "WY")

# Filter rows where the states column contains only CONUS states
huc12_boundaries_conus <- huc12_boundaries[sapply(strsplit(huc12_boundaries$states, ","), function(x) any(x %in% conus_states)), ]

# View filtered results
head(huc12_boundaries_conus)

htdt = readRDS("Variables/huc12_to_downstream_tracts_2021.rds")
htdt_conus <- htdt[htdt$huc12 %in% huc12_boundaries_conus$huc12, ]

saveRDS(htdt_conus, "Variables/htdt_conus.rds")




#-------------------------------------------------------------------------------#
#      Get SA/huc12 lookup table to later mosaic huc12s into service areas      #
#-------------------------------------------------------------------------------#

# Load raster for CRS reference
example_raster <- rast(list.files("Classified HUC12 Rasters", pattern = "\\.tif$", full.names = TRUE)[1])

# Paths to spatial datasets
huc12_shp <- "Variables/huc12_boundaries.gpkg"
sa_shp <- "Service Areas/ServiceAreas_agg.rds"

# Load spatial data
huc12s <- vect(huc12_shp)  
SAs <- readRDS(sa_shp) 

# Ensure CRS consistency
huc12s <- project(huc12s, crs(example_raster))
SAs <- project(SAs, crs(example_raster))

# Define function for intersection. Load everything in the function for parallel processing. 
intersect_sa_huc12 <- function(i) {
  
  library(terra)
  # Paths to spatial datasets
  huc12_shp <- "Variables/huc12_boundaries.gpkg"
  sa_shp <- "Service Areas/ServiceAreas_agg.rds"
  
  # Load spatial data
  huc12s <- vect(huc12_shp)  
  SAs <- readRDS(sa_shp)
  
  
  sa <- SAs[i, ]
  result <- terra::intersect(sa, huc12s)
  if (!is.null(result) && nrow(result) > 0) {
    result_df <- as.data.frame(result)[, c("ID", "huc12")]
    saveRDS(result_df, paste0("Service Areas/Temp lookup/intersection_", i, ".rds"))  # Save intermediate files in this temp folder
    return(result_df)
  }
  return(NULL)
}

plan(multisession, workers = 5)  

# Process in parallel
sa_indices <- seq_len(length(SAs))
results <- future_lapply(sa_indices, intersect_sa_huc12, future.seed = TRUE)

# Combine results from disk instead of memory
all_files <- list.files("Service Areas/Temp lookup", pattern = "intersection_.*\\.rds", full.names = TRUE)
lookup_table <- do.call(rbind, lapply(all_files, readRDS))

# Group by Service Area
lookup_table_grouped <- lookup_table %>%
  group_by(ID) %>%
  summarise(huc12s = list(huc12), .groups = "drop")

# Save the final lookup table
saveRDS(lookup_table_grouped, "Service Areas/Service_Area_HUC12_Lookup.rds")

gc()
print("Process completed with parallel execution and intermediate saving!")

