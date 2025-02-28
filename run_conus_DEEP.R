library(tidyverse)
library(terra)
library(sf)
library(progressr)
library(future.apply)

setwd("C:\\Users\\indumati\\Box\\Paper2_final")


###############

# Create output directory
main_output_dir <- "Classified HUC12 Rasters"
dir.create(main_output_dir, recursive = TRUE, showWarnings = FALSE)

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


# Progress bar
handlers(global = TRUE)  # Enable progress handler globally
handlers("txtprogressbar")  

process_huc12_single = function(huc12_code, output_dir) {
  tryCatch({
    cat(sprintf("\nProcessing unique HUC12: %s\n", huc12_code))
    
    # Load spatial data
    wetland_raster <- rast("Variables/conus_wetlands.tif")
    huc12s_fl <- st_read("Variables/huc12_boundaries.gpkg", quiet = TRUE)
      huc12s_fl <- st_transform(huc12s_fl, st_crs(wetland_raster))
    tracts_centroids <- st_read("Variables/tracts_centroids.shp", quiet = TRUE)
      tracts_centroids <- st_transform(tracts_centroids, st_crs(wetland_raster))
    tracts_data = readRDS("Variables/all_vars_tracts.rds")
      tracts_data <- st_transform(tracts_data, st_crs(wetland_raster))
    huc12_to_downstream_tracts <- readRDS("Variables/htdt_conus.rds")
    
    # Get downstream tracts for this HUC12
    downstream_tracts <- huc12_to_downstream_tracts %>%
      filter(huc12 == huc12_code) %>%
      pull(downstream_tracts) %>%
      unlist()
    
    # Process HUC12 boundary and wetland raster
    huc12_boundary <- huc12s_fl %>%
      filter(huc12 == huc12_code)
    
    wetland_raster_huc12 <- crop(wetland_raster, huc12_boundary) %>%
      mask(huc12_boundary)
    
    
    # Filter centroids and total values
    downstream_centroids <- tracts_centroids %>%
      filter(GEOID %in% downstream_tracts)
    
    tract_data <- tracts_data %>%
      filter(GEOID %in% downstream_tracts) %>%
      arrange(match(GEOID, downstream_tracts))
    
    # Calculate distances from tracts to HUC12 boundary
    distances <- st_distance(downstream_centroids, huc12_boundary) / 1000  # Convert to km
    distances <- as.numeric(distances)
    
    within_40km_indices <- which(distances <= 40)
    if (length(within_40km_indices) == 0) {
      cat("No downstream tracts within 40 km of HUC12\n")
      return(NULL)
    }
    downstream_centroids_filtered <- downstream_centroids[within_40km_indices, ]
    tract_data_filtered <- tract_data[within_40km_indices, ]
    tract_hu_filtered <- tract_data_filtered$housing_units
    tract_hv_filtered <- tract_data_filtered$total_housing_value
    tract_pop_filtered <- tract_data_filtered$pop
    
    tract_hu_filtered <- ifelse(tract_hu_filtered == 0, 0.000001, tract_hu_filtered)
    tract_hv_filtered <- ifelse(tract_hv_filtered == 0, 0.000001, tract_hv_filtered)
    tract_pop_filtered <- ifelse(tract_pop_filtered == 0, 0.000001, tract_pop_filtered)
    
    # Process distances and create decay stacks
    decay_stack <- rast()
    
    for (i in seq_len(nrow(downstream_centroids_filtered))) {
      distance_raster <- distance(wetland_raster_huc12, 
                                  y = downstream_centroids_filtered[i,]) / 1000
      decay_raster <- app(distance_raster, fun = distance_decay)
      decay_stack <- c(decay_stack, decay_raster)
      rm(distance_raster, decay_raster)
      gc()
    }
    
    # Calculate final scores
    score_raster_hu <- terra::weighted.mean(decay_stack, w = tract_hu_filtered, na.rm = TRUE)
    score_raster_hv <- terra::weighted.mean(decay_stack, w = tract_hv_filtered, na.rm = TRUE)
    score_raster_pop <- terra::weighted.mean(decay_stack, w = tract_pop_filtered, na.rm = TRUE)
    
    # Denormalize scores
    unnormalized_score_hu <- score_raster_hu * sum(tract_hu_filtered)
    unnormalized_score_hv <- score_raster_hv * sum(tract_hv_filtered)
    unnormalized_score_pop <- score_raster_pop * sum(tract_pop_filtered)
    
    wetland_raster_huc12[wetland_raster_huc12 == 0] <- NA
    
    # Crop to HUC12 extent
    cropped_hu <- mask(unnormalized_score_hu, wetland_raster_huc12)
    cropped_hv <- mask(unnormalized_score_hv, wetland_raster_huc12)
    cropped_pop <- mask(unnormalized_score_pop, wetland_raster_huc12)
    
    # Save multi-layer raster for this HUC12
    output_file <- file.path(output_dir, sprintf("huc12_%s_raster.tif", huc12_code))
    writeRaster(c(cropped_hu, cropped_hv, cropped_pop), output_file, overwrite = TRUE)
    
    cat(sprintf("\nSaved multi-layer raster for HUC12 %s to: %s\n", huc12_code, output_file))
    
    return(list(huc12_code = huc12_code, file = output_file))
  }, error = function(e) {
    cat(sprintf("Error processing unique HUC12 %s: %s\n", huc12_code, e$message))
    return(NULL)
  })
}



process_batches_startbatch <- function(huc12_file, output_dir, batch_size, start_batch) {
  # Load the HUC12 to downstream tracts data
  huc12_to_downstream_tracts <- readRDS(huc12_file)
  
  unique_huc12s <- unique(huc12_to_downstream_tracts$huc12)
  total_huc12s <- length(unique_huc12s)
  num_batches <- ceiling(total_huc12s / batch_size)
  
  cat(sprintf("Found %d unique HUC12s to process in %d batches (batch size: %d)\n", 
              total_huc12s, num_batches, batch_size))
  
  all_results <- list()  # To collect results across all batches
  
  with_progress({
    for (batch in seq(start_batch, num_batches)) {
      # Stop any old workers before starting a new batch
      future:::ClusterRegistry("stop")  
      plan(multisession, workers = 10)  # Restart workers for new batch
      
      # Determine indices for this batch
      start_idx <- (batch - 1) * batch_size + 1
      end_idx <- min(batch * batch_size, total_huc12s)
      batch_huc12s <- unique_huc12s[start_idx:end_idx]
      
      # Progress message for batch
      cat(sprintf("Processing batch %d of %d (%d HUC12s: %d to %d)\n", 
                  batch, num_batches, length(batch_huc12s), start_idx, end_idx))
      
      # Progress bar for batch
      p <- progressor(along = batch_huc12s)
      
      # Parallel processing for the current batch
      batch_results <- future_lapply(batch_huc12s, function(huc12) {
        p(sprintf("Processing HUC12: %s", huc12))  # Update progress bar
        process_huc12_single(huc12, output_dir)
      }, future.scheduling = 2)  # Improve load balancing
      
      # Filter valid results from the batch and append to all_results
      valid_results <- Filter(Negate(is.null), batch_results)
      all_results <- c(all_results, valid_results)
      
      # Free memory after each batch
      gc()
    }
  })
  
  return(all_results)
}


output_dir <- "Classified HUC12 Rasters"
huc12_file <- "Variables/htdt_conus.rds"

plan(multisession, workers=10)
# Started 9:00am 2/20
results <- process_batches_startbatch(
  huc12_file, 
  output_dir, 
  batch_size = 250, 
  start_batch = 1
)


######## If it crashes, exclude already processed files #######


output_dir <- "Classified HUC12 Rasters"
huc12_file <- "Variables/htdt_conus.rds"

# List already processed HUC12 files
processed_files <- list.files("Classified HUC12 Rasters", pattern = "huc12_.*_raster\\.tif$", full.names = TRUE)

# Extract HUC12 codes from filenames
processed_huc12s <- gsub("Classified HUC12 Rasters/huc12_|_raster\\.tif", "", processed_files)

# Load the original data frame
huc12_to_downstream_tracts <- readRDS("Variables/htdt_conus.rds")

# Filter out already processed HUC12s (keep only rows with unprocessed HUC12s)
remaining_huc12s <- huc12_to_downstream_tracts %>%
  filter(!(huc12 %in% processed_huc12s))

saveRDS(remaining_huc12s, "Variables/remaining.rds")
plan(multisession, workers=10)
# Run sequentially
results <- process_batches_startbatch(
  huc12_file = "Variables/remaining.rds",  # Use filtered HUC12s
  output_dir = "Classified HUC12 Rasters", 
  batch_size = 250, 
  start_batch = 1
)

