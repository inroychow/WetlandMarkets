library(terra)
library(stringr)
library(tidyverse)
library(terra)
library(sf)
library(progressr)
library(future.apply)

setwd("C:/Users/indumati/Box/Paper2_final")

############################################################################################################################
#############                                                                                                     ##########
############# IDENTIFY THE HUC12s THAT ARE MISSING, AND GET THE WETLAND RASTER FOR THEM. SAVE THESE IN THE FOLDER ##########
#############                           "Classified HUC12 Rasters/Missing HUC12s - 40km"                          ##########
#############                                                                                                     ##########
############################################################################################################################


# 1. Load HUC12-to-Downstream Tracts Data
htdt_conus <- readRDS("Variables/htdt_conus.rds")

# Extract all unique HUC12s that should be processed
expected_huc12s <- unique(htdt_conus$huc12)

# 2. Identify which HUC12 rasters already exist
huc12_rasters_folder <- "Classified HUC12 Rasters"
processed_files <- list.files(huc12_rasters_folder, 
                              pattern="^huc12_.*_raster\\.tif$", 
                              full.names=FALSE)

# Extract HUC12 codes from the filenames
processed_huc12s <- str_extract(processed_files, "(?<=huc12_)[0-9]+(?=_raster\\.tif)")

# 3. Find missing HUC12s
missing_huc12s <- setdiff(expected_huc12s, processed_huc12s)

# 4. Load the HUC12 boundaries and the main wetland raster
huc12_boundaries <- vect("Variables/huc12_boundaries.gpkg")
wetland_raster   <- rast("Variables/conus_wetlands.tif")

# 5. Create output folder for extracted wetlands
extracted_wetlands_folder <- "Classified HUC12 Rasters/Missing HUC12 Rasters - 40km"
dir.create(extracted_wetlands_folder, recursive = TRUE, showWarnings = FALSE)

# 6. Extract wetland rasters for each missing HUC12
for (h12 in missing_huc12s) {
  cat("Extracting wetlands for missing HUC12:", h12, "\n")
  
  # --- Subset the SpatVector directly
  h12_boundary <- huc12_boundaries[huc12_boundaries$huc12 == h12, ]
  
  # If nothing is returned, skip
  if (nrow(h12_boundary) == 0) {
    cat("huc12", h12, "not found in boundary file, skipping.\n")
    next
  }
  
  # Crop and mask the wetlands raster
  wetland_h12 <- crop(wetland_raster, h12_boundary)
  wetland_h12 <- mask(wetland_h12, h12_boundary)
  
  # Save the wetland raster
  out_file <- file.path(extracted_wetlands_folder, 
                        paste0("huc12_", h12, "_wetlands.tif"))
  writeRaster(wetland_h12, out_file, overwrite=TRUE)
  
  cat("Saved wetland raster for missing HUC12:", h12, "\n")
}

cat("Done extracting wetlands for all missing HUC12s.\n")



#######################################################################################################################
#############                                                                                                ##########
############# Run Distance Decay Classification workflow on the "missing" huc12s based on the tract BOUNDARY ##########
#############                              rather than the tract centroid.                                   ##########
#############                                                                                                ##########
########################################################################################################################


###############

setwd("L:/Wetland Flood Mitigation")

handlers(global = TRUE)  # Enable progress handler globally
handlers("txtprogressbar")  



# Func to sanitize service area names for file paths
sanitize_path <- function(path) {
  path <- gsub("[[:space:]]+", "_", path)
  path <- gsub("[^[:alnum:]_/-]", "", path)
  return(path)
}

# Distance decay function
distance_decay <- function(distance_km) {
  distance_km_values <- c(0, 10, 20, 30, 40)
  decay_values <- c(0.25, 1, 0.75, 1, 0)
  
  interp_values <- approx(
    x = distance_km_values,
    y = decay_values,
    xout = distance_km,
    rule = 2,
    yleft = 0,
    yright = 0,
    method = "linear"
  )$y
  return(interp_values)
}

process_huc12_single <- function(huc12_code, output_dir) {
  tryCatch({
    cat(sprintf("\nProcessing HUC12: %s\n", huc12_code))
    
    # Load spatial data
    wetland_raster <- rast("Variables/conus_wetlands.tif")
    huc12s_fl <- st_read("Variables/huc12_boundaries.gpkg", quiet = TRUE)
    huc12s_fl <- st_transform(huc12s_fl, st_crs(wetland_raster))  # Fixed missing )
    
    tracts_sf <- readRDS("Variables/all_vars_tracts.rds")  # Tract boundaries
    tracts_sf <- st_transform(tracts_sf, st_crs(wetland_raster))  # Fixed missing )
    
    huc12_to_downstream_tracts <- readRDS("Variables/htdt_conus.rds")
    
    # Get downstream tracts for this HUC12
    downstream_tracts <- huc12_to_downstream_tracts %>%
      filter(huc12 == huc12_code) %>%
      pull(downstream_tracts) %>%
      unlist()
    
    # Get HUC12 boundary
    huc12_boundary <- huc12s_fl %>%
      filter(huc12 == huc12_code)
    
    # Get downstream tract polygons
    downstream_tracts_sf <- tracts_sf %>%
      filter(GEOID %in% downstream_tracts)
    
    # Compute distances from tracts to HUC12 boundary
    distances <- st_distance(downstream_tracts_sf, huc12_boundary) / 1000  # Convert to km
    distances <- apply(distances, 1, min)  # Get minimum distance for each tract
    
    # Filter tracts within 40 km
    within_40km_indices <- which(distances <= 40)
    if (length(within_40km_indices) == 0) {
      cat("No downstream tracts within 40 km\n")
      return(NULL)
    }
    
    downstream_tracts_filtered <- downstream_tracts_sf[within_40km_indices, ]
    
    # Load housing data
    tract_data <- readRDS("Variables/all_vars_tracts.rds") %>%
      filter(GEOID %in% downstream_tracts_filtered$GEOID)
    
    tract_hu <- tract_data$housing_units
    tract_hv <- tract_data$property_value
    tract_pop <- tract_data$pop
    
    # Avoid division by zero
    tract_hu <- ifelse(tract_hu == 0, 0.000001, tract_hu)
    tract_hv <- ifelse(tract_hv == 0, 0.000001, tract_hv)
    tract_pop <- ifelse(tract_pop == 0, 0.000001, tract_pop)
    
    # Crop wetland raster to HUC12
    wetland_raster_huc12 <- crop(wetland_raster, huc12_boundary) %>%
      mask(huc12_boundary)
    
    # Create decay raster stack
    decay_stack <- rast()
    
    for (i in seq_len(nrow(downstream_tracts_filtered))) {
      distance_raster <- distance(wetland_raster_huc12, 
                                  y = downstream_tracts_filtered[i,]) / 1000
      decay_raster <- app(distance_raster, fun = distance_decay)
      decay_stack <- c(decay_stack, decay_raster)
      rm(distance_raster, decay_raster)
      gc()
    }
    
    # Compute weighted decay
    score_raster_hu <- terra::weighted.mean(decay_stack, w = tract_hu, na.rm = TRUE)
    score_raster_hv <- terra::weighted.mean(decay_stack, w = tract_hv, na.rm = TRUE)
    score_raster_pop <- terra::weighted.mean(decay_stack, w = tract_pop, na.rm = TRUE)
    
    # Denormalize scores
    unnormalized_score_hu <- score_raster_hu * sum(tract_hu)
    unnormalized_score_hv <- score_raster_hv * sum(tract_hv)
    unnormalized_score_pop <- score_raster_pop * sum(tract_pop)
    
    # Mask final rasters to HUC12
    cropped_hu <- mask(unnormalized_score_hu, wetland_raster_huc12)
    cropped_hv <- mask(unnormalized_score_hv, wetland_raster_huc12)
    cropped_pop <- mask(unnormalized_score_pop, wetland_raster_huc12)
    
    # Save multi-layer raster
    output_file <- file.path(output_dir, sprintf("huc12_%s_raster.tif", huc12_code))
    writeRaster(c(cropped_hu, cropped_hv, cropped_pop), output_file, overwrite = TRUE)
    
    cat(sprintf("\nSaved raster for HUC12 %s: %s\n", huc12_code, output_file))
    
    return(list(huc12_code = huc12_code, file = output_file))
  }, error = function(e) {
    cat(sprintf("Error processing HUC12 %s: %s\n", huc12_code, e$message))
    return(NULL)
  })
}




plan(multisession, workers=10)
# Started 9:00am 2/20
results <- process_batches_startbatch(
  huc12_file = "Variables/remaining.rds", 
  output_dir = "Classified HUC12 Rasters/40km_tract_boundary_classified", 
  batch_size = 1, 
  start_batch = 1
)

#  030502080609 is the only huc12 that is not running.


####################################################################################
#######     Impute zeros for the huc12s that STILL have no downstream tracts  ######
####################################################################################
library(stringr)

# Read the list of remaining HUC12s that were originally missing
remaining_huc12s <- readRDS("Variables/remaining.rds")
remaining_huc12s = remaining_huc12s$huc12
# Get HUC12s that were successfully classified (from filenames)
classified_folder <- "Classified HUC12 Rasters/40km_tract_boundary_classified"
classified_files <- list.files(classified_folder, pattern="^huc12_.*_raster\\.tif$", full.names=FALSE)

# Extract HUC12 codes from filenames
classified_huc12s <- str_extract(classified_files, "(?<=huc12_)[0-9]+(?=_raster\\.tif)")

# Subtract the classified HUC12s from the remaining ones
still_missing_huc12s <- setdiff(remaining_huc12s, classified_huc12s)

# Save the list of still-missing HUC12s
saveRDS(still_missing_huc12s, "Variables/still_missing_huc12s.rds")

cat("Number of HUC12s still missing:", length(still_missing_huc12s), "\n")


####################

library(terra)

# Load the list of still-missing HUC12s
still_missing_huc12s <- readRDS("Variables/still_missing_huc12s.rds")

# Define folders
missing_folder <- "Classified HUC12 Rasters/Missing HUC12s - 40km"
classified_folder <- "Classified HUC12 Rasters/40km_tract_boundary_classified"

# Loop over each missing HUC12
for (h12 in still_missing_huc12s) {
  cat("Processing HUC12:", h12, "\n")
  
  # Path to input wetland raster
  wetland_file <- file.path(missing_folder, paste0("huc12_", h12, "_wetlands.tif"))
  
  if (file.exists(wetland_file)) {
    # Load raster
    wetland_raster <- rast(wetland_file)
    
    # Create a zero raster, but preserve NA values
    zero_raster <- wetland_raster
    zero_raster[!is.na(wetland_raster)] <- 0  # Assign 0 to wetland cells only
    
    # Create a 3-layer raster (for HU, HV, and Pop) while keeping NA values
    zero_stack <- c(zero_raster, zero_raster, zero_raster)
    
    # Define output path
    output_file <- file.path(classified_folder, paste0("huc12_", h12, "_raster.tif"))
    
    # Save the modified raster
    writeRaster(zero_stack, output_file, overwrite=TRUE, datatype="FLT4S")
    
    cat("Saved zero-assigned raster for wetlands in HUC12:", h12, "\n")
  } else {
    cat("Wetland raster missing for HUC12:", h12, " - Skipping\n")
  }
}

cat("\n*** Zero raster assignment complete for missing HUC12s ***\n")
