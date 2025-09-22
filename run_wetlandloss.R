# ---------------------------------------------------------------
# ---------------------------------------------------------------
#  
#
#           Classify wetland loss cells (1985–2021) by upstream
#           flood protection value using a distance decay function
#           and tract-level weights (housing units, value, pop).
# ---------------------------------------------------------------

library(tidyverse)
library(terra)
library(sf)
library(progressr)
library(future.apply)

# --------------------- CONFIGURE ROOT PATH ---------------------

datadir <- Sys.getenv("DATA_DIR", "C:/Users/indumati/Box/Paper2_final")
output_dir <- file.path(datadir, "Wetland Loss", "Classified HUC12 Rasters - Loss")
error_log <- file.path(datadir, "Wetland Loss", "huc12_error_log.rds")
no_centroids_log <- file.path(datadir, "Wetland Loss", "huc12_no_centroids_within_40km.rds")

# Ensure required output/log folders exist
dir.create(output_dir, recursive = TRUE, showWarnings = FALSE)

# --------------------- UTILITY FUNCTIONS ------------------------

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
}

handlers(global = TRUE)
handlers("txtprogressbar")

# --------------------- Single HUC12 process function -------------------

process_huc12_single <- function(huc12_code, output_dir,
                                 error_log = error_log,
                                 no_centroids_log = no_centroids_log) {
  tryCatch({
    cat(sprintf("\nProcessing unique HUC12: %s\n", huc12_code))
    
    wetland_loss_raster <- rast(file.path(datadir, "Wetland Loss", "wetland_to_dev.tif"))
    huc12s <- st_read(file.path(datadir, "Variables", "huc12_boundaries.gpkg"), quiet = TRUE) %>%
      st_transform(crs(wetland_loss_raster))
    tracts_centroids <- st_read(file.path(datadir, "Variables", "tracts_centroids.shp"), quiet = TRUE) %>%
      st_transform(crs(wetland_loss_raster))
    tracts_data <- readRDS(file.path(datadir, "Variables", "all_vars_tracts.rds")) %>%
      st_transform(crs(wetland_loss_raster))
    huc12_to_downstream_tracts <- readRDS(file.path(datadir, "Variables", "htdt_conus.rds"))
    
    downstream_tracts <- huc12_to_downstream_tracts %>%
      filter(huc12 == huc12_code) %>%
      pull(downstream_tracts) %>%
      unlist()
    
    huc12_boundary <- huc12s %>% filter(huc12 == huc12_code)
    lost_wetlands_huc12 <- crop(wetland_loss_raster, huc12_boundary) %>% mask(huc12_boundary)
    
    if (all(is.na(values(lost_wetlands_huc12)))) {
      cat(sprintf("No lost wetlands in HUC12 %s, skipping...\n", huc12_code))
      return(NULL)
    }
    
    downstream_centroids <- tracts_centroids %>% filter(GEOID %in% downstream_tracts)
    tract_data <- tracts_data %>% filter(GEOID %in% downstream_tracts) %>% arrange(match(GEOID, downstream_tracts))
    
    distances <- as.numeric(st_distance(downstream_centroids, huc12_boundary) / 1000)
    within_40km <- which(distances <= 40)
    if (length(within_40km) == 0) {
      cat(sprintf("No downstream tracts within 40 km for HUC12 %s\n", huc12_code))
      log_list <- if (file.exists(no_centroids_log)) readRDS(no_centroids_log) else list()
      log_list[[huc12_code]] <- list(huc12 = huc12_code, reason = "No centroids within 40km")
      saveRDS(log_list, no_centroids_log)
      return(NULL)
    }
    
    downstream_centroids <- downstream_centroids[within_40km, ]
    tract_data <- tract_data[within_40km, ]
    tract_hu <- ifelse(tract_data$housing_units == 0, 1e-6, tract_data$housing_units)
    tract_hv <- ifelse(tract_data$total_housing_value == 0, 1e-6, tract_data$total_housing_value)
    tract_pop <- ifelse(tract_data$pop == 0, 1e-6, tract_data$pop)
    
    decay_stack <- rast()
    for (i in seq_len(nrow(downstream_centroids))) {
      d_raster <- distance(lost_wetlands_huc12, downstream_centroids[i,]) / 1000
      decay_raster <- app(d_raster, fun = distance_decay)
      decay_stack <- c(decay_stack, decay_raster)
      rm(d_raster, decay_raster); gc()
    }
    
    score_hu <- weighted.mean(decay_stack, w = tract_hu, na.rm = TRUE)
    score_hv <- weighted.mean(decay_stack, w = tract_hv, na.rm = TRUE)
    score_pop <- weighted.mean(decay_stack, w = tract_pop, na.rm = TRUE)
    
    lost_wetlands_huc12[lost_wetlands_huc12 == 0] <- NA
    cropped_hu <- mask(score_hu * sum(tract_hu), lost_wetlands_huc12)
    cropped_hv <- mask(score_hv * sum(tract_hv), lost_wetlands_huc12)
    cropped_pop <- mask(score_pop * sum(tract_pop), lost_wetlands_huc12)
    
    out_file <- file.path(output_dir, sprintf("huc12_%s_raster1985loss.tif", huc12_code))
    writeRaster(c(cropped_hu, cropped_hv, cropped_pop), out_file, overwrite = TRUE)
    cat(sprintf("\nSaved raster for HUC12 %s to: %s\n", huc12_code, out_file))
    return(list(huc12_code = huc12_code, file = out_file))
  }, error = function(e) {
    cat(sprintf("Error processing HUC12 %s: %s\n", huc12_code, e$message))
    log_list <- if (file.exists(error_log)) readRDS(error_log) else list()
    log_list[[huc12_code]] <- e$message
    saveRDS(log_list, error_log)
    return(NULL)
  })
}

# --------------------- Batch process function --------------------------

process_batches_startbatch <- function(huc12_file, output_dir, batch_size, start_batch,
                                       log_file = file.path(datadir, "Wetland Loss", "batch_progress_log.txt")) {
  htdt <- readRDS(huc12_file)
  unique_huc12s <- unique(htdt$huc12)
  total <- length(unique_huc12s)
  num_batches <- ceiling(total / batch_size)
  
  cat(sprintf("Found %d HUC12s in %d batches\n", total, num_batches))
  log_conn <- file(log_file, open = "a")
  all_results <- list()
  
  with_progress({
    for (batch in seq(start_batch, num_batches)) {
      future:::ClusterRegistry("stop")
      plan(multisession, workers = 10)
      start_idx <- (batch - 1) * batch_size + 1
      end_idx <- min(batch * batch_size, total)
      batch_huc12s <- unique_huc12s[start_idx:end_idx]
      
      msg <- sprintf("Batch %d of %d (%d HUC12s) at %s\n", batch, num_batches,
                     length(batch_huc12s), Sys.time())
      cat(msg); writeLines(msg, log_conn)
      
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
  close(log_conn)
  return(all_results)
}

# --------------------- Run main batch ---------------------------

plan(multisession, workers = 10)
results <- process_batches_startbatch(
  huc12_file = file.path(datadir, "Variables", "htdt_conus.rds"),
  output_dir = output_dir,
  batch_size = 250,
  start_batch = 1
)

# ------------------ filter remaining and run again if needed ------------------

processed <- list.files(output_dir, pattern = "huc12_.*_raster1985loss\\.tif$", full.names = TRUE)
processed_ids <- gsub(".*huc12_(.*)_raster1985loss\\.tif$", "\\1", processed)
huc12_df <- readRDS(file.path(datadir, "Variables", "htdt_conus.rds"))
remaining <- huc12_df %>% filter(!(huc12 %in% processed_ids))
saveRDS(remaining, file.path(datadir, "Variables", "remaining.rds"))

plan(multisession, workers = 10)
results <- process_batches_startbatch(
  huc12_file = file.path(datadir, "Variables", "remaining.rds"),
  output_dir = file.path(datadir, "Wetland Loss", "Output loss - 3.5"),
  batch_size = 500,
  start_batch = 1
)
