
#           Classify wetland cells by upstream flood protection value
#           using a distance decay function and tract-level weights
#           (housing units, value, population).
# ---------------------------------------------------------------

library(tidyverse)
library(terra)
library(sf)
library(progressr)
library(future.apply)

# --------------------- Set root path (change as needed) ---------
Sys.setenv("DATA_DIR" = "C:/Users/indumati/Box/Paper2_final")

datadir <- Sys.getenv("DATA_DIR", "C:/Users/indumati/Box/Paper2_final")
output_dir <- file.path(datadir, "Classified HUC12 Rasters")

# Ensure output folder exists
dir.create(output_dir, recursive = TRUE, showWarnings = FALSE)

# --------------------- Utility functions ------------------------

sanitize_path <- function(path) {
  path <- gsub("[[:space:]]+", "_", path)
  path <- gsub("[^[:alnum:]_/-]", "", path)
  return(path)
}

distance_decay <- function(distance_km) {
  distance_km_values <- c(0, 10, 20, 30, 40)
  decay_values <- c(0.25, 1, 0.75, 1, 0)
  approx(x = distance_km_values, y = decay_values, xout = distance_km,
         rule = 2, yleft = 0, yright = 0, method = "linear")$y
} #Interpolate values between

handlers(global = TRUE)
handlers("txtprogressbar")

# --------------------- Processing function for one huc12 ----------------------

process_huc12_single <- function(huc12_code, output_dir) {
  tryCatch({
    cat(sprintf("\nProcessing unique HUC12: %s\n", huc12_code))
    
    wetland_raster <- rast(file.path(datadir, "Variables", "conus_wetlands.tif"))
    huc12s_fl <- st_read(file.path(datadir, "Variables", "huc12_boundaries.gpkg"), quiet = TRUE) %>%
      st_transform(crs(wetland_raster))
    tracts_centroids <- st_read(file.path(datadir, "Variables", "tracts_centroids.shp"), quiet = TRUE) %>%
      st_transform(crs(wetland_raster))
    tracts_data <- readRDS(file.path(datadir, "Variables", "all_vars_tracts.rds")) %>%
      st_transform(crs(wetland_raster))
    huc12_to_downstream_tracts <- readRDS(file.path(datadir, "Variables", "htdt_conus.rds"))
    
    downstream_tracts <- huc12_to_downstream_tracts %>%
      filter(huc12 == huc12_code) %>%
      pull(downstream_tracts) %>%
      unlist()
    
    huc12_boundary <- huc12s_fl %>% filter(huc12 == huc12_code)
    wetland_raster_huc12 <- crop(wetland_raster, huc12_boundary) %>% mask(huc12_boundary)
    
    downstream_centroids <- tracts_centroids %>% filter(GEOID %in% downstream_tracts)
    tract_data <- tracts_data %>% filter(GEOID %in% downstream_tracts) %>% arrange(match(GEOID, downstream_tracts))
    
    distances <- as.numeric(st_distance(downstream_centroids, huc12_boundary) / 1000)
    within_40km_indices <- which(distances <= 40)
    if (length(within_40km_indices) == 0) {
      cat("No downstream tracts within 40 km of HUC12\n")
      return(NULL)
    }
    
    downstream_centroids_filtered <- downstream_centroids[within_40km_indices, ]
    tract_data_filtered <- tract_data[within_40km_indices, ]
    tract_hu <- ifelse(tract_data_filtered$housing_units == 0, 1e-6, tract_data_filtered$housing_units)
    tract_hv <- ifelse(tract_data_filtered$total_housing_value == 0, 1e-6, tract_data_filtered$total_housing_value)
    tract_pop <- ifelse(tract_data_filtered$pop == 0, 1e-6, tract_data_filtered$pop)
    
    decay_stack <- rast()
    for (i in seq_len(nrow(downstream_centroids_filtered))) {
      distance_raster <- distance(wetland_raster_huc12, downstream_centroids_filtered[i,]) / 1000
      decay_raster <- app(distance_raster, fun = distance_decay)
      decay_stack <- c(decay_stack, decay_raster)
      rm(distance_raster, decay_raster); gc()
    }
    
    score_raster_hu <- weighted.mean(decay_stack, w = tract_hu, na.rm = TRUE)
    score_raster_hv <- weighted.mean(decay_stack, w = tract_hv, na.rm = TRUE)
    score_raster_pop <- weighted.mean(decay_stack, w = tract_pop, na.rm = TRUE)
    
    wetland_raster_huc12[wetland_raster_huc12 == 0] <- NA
    cropped_hu <- mask(score_raster_hu * sum(tract_hu), wetland_raster_huc12)
    cropped_hv <- mask(score_raster_hv * sum(tract_hv), wetland_raster_huc12)
    cropped_pop <- mask(score_raster_pop * sum(tract_pop), wetland_raster_huc12)
    
    output_file <- file.path(output_dir, sprintf("huc12_%s_raster.tif", huc12_code))
    writeRaster(c(cropped_hu, cropped_hv, cropped_pop), output_file, overwrite = TRUE)
    cat(sprintf("\nSaved raster for HUC12 %s to: %s\n", huc12_code, output_file))
    return(list(huc12_code = huc12_code, file = output_file))
  }, error = function(e) {
    cat(sprintf("Error processing HUC12 %s: %s\n", huc12_code, e$message))
    return(NULL)
  })
}

# --------------------- Batch wrapper for parallel execution  ----------------------------

process_batches_startbatch <- function(huc12_file, output_dir, batch_size, start_batch) {
  htdt <- readRDS(huc12_file)
  unique_huc12s <- unique(htdt$huc12)
  total_huc12s <- length(unique_huc12s)
  num_batches <- ceiling(total_huc12s / batch_size)
  
  cat(sprintf("Found %d unique HUC12s in %d batches (batch size: %d)\n", 
              total_huc12s, num_batches, batch_size))
  
  all_results <- list()
  with_progress({
    for (batch in seq(start_batch, num_batches)) {
      future:::ClusterRegistry("stop")
      plan(multisession, workers = 10)
      
      start_idx <- (batch - 1) * batch_size + 1
      end_idx <- min(batch * batch_size, total_huc12s)
      batch_huc12s <- unique_huc12s[start_idx:end_idx]
      
      cat(sprintf("Processing batch %d of %d (%d HUC12s: %d to %d)\n", 
                  batch, num_batches, length(batch_huc12s), start_idx, end_idx))
      
      p <- progressor(along = batch_huc12s)
      batch_results <- future_lapply(batch_huc12s, function(huc12) {
        p(sprintf("Processing HUC12: %s", huc12))
        process_huc12_single(huc12, output_dir)
      }, future.scheduling = 2)
      
      valid_results <- Filter(Negate(is.null), batch_results)
      all_results <- c(all_results, valid_results)
      gc()
    }
  })
  return(all_results)
}

# --------------------- Run  -------------------------------

plan(multisession, workers = 10)
results <- process_batches_startbatch(
  huc12_file = file.path(datadir, "Variables", "htdt_conus.rds"),
  output_dir = output_dir,
  batch_size = 250,
  start_batch = 1
)

# ------- If it crashes, Run again by filtering remaining huc12s -------

processed_files <- list.files(output_dir, pattern = "huc12_.*_raster\\.tif$", full.names = TRUE)
processed_huc12s <- gsub(".*huc12_(.*)_raster\\.tif$", "\\1", processed_files)
htdt <- readRDS(file.path(datadir, "Variables", "htdt_conus.rds"))
remaining_huc12s <- htdt %>% filter(!(huc12 %in% processed_huc12s))
saveRDS(remaining_huc12s, file.path(datadir, "Variables", "remaining.rds"))

plan(multisession, workers = 10)
results <- process_batches_startbatch(
  huc12_file = file.path(datadir, "Variables", "remaining.rds"),
  output_dir = output_dir,
  batch_size = 250,
  start_batch = 1
)
