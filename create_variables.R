
data_dir <- Sys.getenv("DATA_DIR", "C:/path/to/project")

# ---- 2. Sub‑directories ----
geo_dir         <- file.path(data_dir, "Geo Boundaries")
raster_dir      <- file.path(data_dir, "LCMAP Rasters")
service_area_dir<- file.path(data_dir, "Service Areas")
output_dir      <- file.path(data_dir, "Variables", "State-by-state")
classified_dir  <- file.path(data_dir, "Classified HUC12 Rasters")
python_dir      <- file.path(data_dir, "Python")
temp_lookup_dir <- file.path(service_area_dir, "Temp lookup")

# Create persistent folders if they do not exist
invisible(lapply(c(output_dir, temp_lookup_dir), dir.create,
                 recursive = TRUE, showWarnings = FALSE))

# -------------------- LIBRARIES --------------------------------
suppressPackageStartupMessages({
  library(tidyverse)
  library(data.table)
  library(xml2)
  library(terra)
  library(sf)
  library(tmap)
  library(tidycensus)
  library(future.apply)   # for parallel loops
})

# -------------------- INPUT DATA -------------------------------

# Wetland raster (2021 LCMAP primary land‑cover)
wetland_raster <- rast(file.path(raster_dir,
                                 "LCMAP_CU_2021_V13_LCPRI.tif"))

# HUC12 boundaries 
huc12_boundaries <- st_read(file.path(geo_dir, "HUC12s",
                                      "WBD_National_GDB.gdb"),
                            layer = "WBDHU12", quiet = TRUE) %>%
  st_make_valid() %>%
  st_transform(crs(wetland_raster))

# CONUS state boundaries (as SpatVector for terra)
conus <- readRDS(file.path(geo_dir, "USA States", "conus.rds")) %>%
  project(crs(wetland_raster))

# Service‑area polygons (SpatVector)
serviceareas <- vect(file.path(service_area_dir,
                               "ServiceAreasShapefile.shp")) %>%
  project(crs(wetland_raster))

# Upstream HUC12 → downstream tract lookup from Python (data.table)
upstream_df <- fread(file.path(python_dir,
                               "service_area_huc12_to_downstream_tracts.csv"))

# -------------------- STATE LOOP -------------------------------

states <- unique(conus$STATEFP)

for (state in states) {
  message("\n── Processing state: ", state)
  
  ## --- Spatial subsets -------------------------------------------------
  state_boundary <- conus[conus$STATEFP == state, ]
  state_SAs      <- mask(serviceareas, state_boundary)
  state_SAs$service_area <- gsub("^Document/", "", state_SAs$FolderPath)
  
  # Crop wetland raster
  state_wetlands <- crop(wetland_raster, state_boundary)
  
  ## --- Census tract attributes ----------------------------------------
  get_acs_tract <- function(var){
    get_acs(geography = "tract", variables = var, state = state,
            year = 2021, geometry = TRUE) %>%
      st_transform(crs(state_wetlands))
  }
  hv_tracts  <- get_acs_tract("B25077_001") %>% rename(property_value = estimate)
  hu_tracts  <- get_acs_tract("B25001_001") %>% rename(housing_units  = estimate)
  pop_tracts <- get_acs_tract("B01003_001") %>% rename(population     = estimate)
  
  ## --- HUC‑12 raster (per state) --------------------------------------
  state_huc12_boundaries <- st_crop(huc12_boundaries, state_boundary) %>%
    st_transform(crs(state_wetlands))
  
  state_huc12_vect   <- vect(state_huc12_boundaries)
  state_huc12_raster <- rasterize(state_huc12_vect, state_wetlands,
                                  field = "huc12", touches = TRUE)
  state_wetland_huc12_raster <- mask(state_huc12_raster, state_wetlands)
  
  ## --- Tract centroids -------------------------------------------------
  tracts_centroids <- st_centroid(hu_tracts) %>%
    select(GEOID, housing_units) %>% st_as_sf()
  
  ## --- Upstream → downstream look‑up (state subset) -------------------
  upstream_df <- upstream_df %>% mutate(huc12 = str_pad(huc12, 12, "0"))
  
  upstream_subset <- upstream_df %>%
    filter(huc12 %in% state_huc12_boundaries$huc12)
  
  state_htdt <- upstream_subset %>%
    mutate(downstream_tracts = str_remove_all(downstream_tracts, "[\\[\\]\']") %>%
             str_split(",\\s*")) %>%
    group_by(huc12) %>%
    summarise(downstream_tracts = list(unique(unlist(downstream_tracts))),
              .groups = "drop") %>%
    mutate(huc12 = str_pad(huc12, 12, "0"),
           downstream_tracts = lapply(downstream_tracts,
                                      function(x) str_pad(x, 11, "0")))
  
  ## --- Per‑state output -----------------------------------------------
  state_dir <- file.path(output_dir, state)
  dir.create(state_dir, recursive = TRUE, showWarnings = FALSE)
  
  writeRaster(state_wetlands,            file.path(state_dir, "wetlands.tif"),  overwrite = TRUE)
  writeRaster(state_huc12_raster,        file.path(state_dir, "huc12_raster.tif"), overwrite = TRUE)
  st_write(tracts_centroids,             file.path(state_dir, "tracts_centroids.shp"), delete_dsn = TRUE, quiet = TRUE)
  st_write(state_huc12_boundaries,       file.path(state_dir, "huc12_boundaries.gpkg"), append = FALSE, quiet = TRUE)
  saveRDS(state_htdt,                    file.path(state_dir, "huc12_to_downstream_tracts.rds"))
  saveRDS(hu_tracts,                     file.path(state_dir, "hu_tracts"))
  saveRDS(hv_tracts,                     file.path(state_dir, "hv_tracts"))
  saveRDS(pop_tracts,                    file.path(state_dir, "pop_tracts"))
  
  message("✔ Completed: ", state)
}

# -------------------- COMBINE STATE FILES ------------------------------

message("\n── Combining per‑state outputs (may take time)…")

combine_state_rds <- function(fname) {
  do.call(bind_rows, lapply(states, function(st) {
    readRDS(file.path(output_dir, st, fname))
  }))
}

# huc12 → downstream tracts (all states)
all_htdt <- combine_state_rds("huc12_to_downstream_tracts.rds") %>%
  distinct()

# Tract centroids (sf)
all_tract_centroids <- do.call(bind_rows, lapply(states, function(st) {
  st_read(file.path(output_dir, st, "tracts_centroids.shp"), quiet = TRUE)
}))

# Attributes
all_hu  <- combine_state_rds("hu_tracts")
all_hv  <- combine_state_rds("hv_tracts")
all_pop <- combine_state_rds("pop_tracts")

# -------------------- HUC12 FILTER (CONUS) -----------------------------

huc12_boundaries <- st_read(file.path(data_dir, "Variables", "huc12_boundaries.gpkg"), quiet = TRUE)

conus_states <- c("AL","AZ","AR","CA","CO","CT","DE","FL","GA","ID","IL","IN","IA","KS","KY","LA","ME","MD","MA","MI","MN","MS","MO","MT","NE","NV","NH","NJ","NM","NY","NC","ND","OH","OK","OR","PA","RI","SC","SD","TN","TX","UT","VT","VA","WA","WV","WI","WY")

huc12_boundaries_conus <- huc12_boundaries[sapply(strsplit(huc12_boundaries$states, ","),
                                                  function(x) any(x %in% conus_states)), ]

htdt_conus <- all_htdt[all_htdt$huc12 %in% huc12_boundaries_conus$huc12, ]

saveRDS(htdt_conus, file.path(data_dir, "Variables", "htdt_conus.rds"))

# -------------------- SA → HUC12 LOOK‑UP -------------------------------

example_raster <- rast(list.files(classified_dir, pattern = "\\.tif$", full.names = TRUE)[1])

huc12_shp <- file.path(data_dir, "Variables", "huc12_boundaries.gpkg")
sa_rds    <- file.path(service_area_dir, "ServiceAreas_agg.rds")

huc12s <- vect(huc12_shp)     %>% project(crs(example_raster))
SAs    <- readRDS(sa_rds)      %>% project(crs(example_raster))

intersect_sa_huc12 <- function(i) {
  library(terra)
  huc12s <- vect(huc12_shp)
  SAs    <- readRDS(sa_rds)
  sa     <- SAs[i, ]
  res    <- terra::intersect(sa, huc12s)
  if (!is.null(res) && nrow(res) > 0) {
    res_df <- as.data.frame(res)[, c("ID", "huc12")]
    saveRDS(res_df, file.path(temp_lookup_dir, sprintf("intersection_%s.rds", i)))
    return(res_df)
  }
  NULL
}

plan(multisession, workers = 5)
sa_indices <- seq_len(length(SAs))
future_lapply(sa_indices, intersect_sa_huc12, future.seed = TRUE)

# Merge intermediate files → final lookup
lookup_table <- do.call(rbind, lapply(list.files(temp_lookup_dir,
                                                 pattern = "intersection_.*\\.rds",
                                                 full.names = TRUE),
                                      readRDS)) %>%
  group_by(ID) %>%
  summarise(huc12s = list(huc12), .groups = "drop")

saveRDS(lookup_table, file.path(service_area_dir, "Service_Area_HUC12_Lookup.rds"))

message("\n Process completed.")
