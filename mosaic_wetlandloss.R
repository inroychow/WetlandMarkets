library(terra)
library(sf)
library(dplyr)
library(future.apply)
library(fst)

setwd("C:/Users/indumati/Box/Paper2_final")

#--- Define workers (adjust as needed)
plan(multisession, workers = 5)  

# Paths
lookup_table_path <- "Service Areas/Service_Area_HUC12_Lookup.rds"
output_dir        <- "Service Area Mosaics Loss"
raster_dir        <- "Wetland Loss/Output loss 2021"
sa_polygons       <- readRDS("Service Areas/ServiceAreas_agg.rds")

# Load the lookup table
lookup_table <- readRDS(lookup_table_path)

# Subset only to service areas that have at least one HUC12 with wetland loss
huc12_files <- list.files(raster_dir, pattern = "tif$", full.names = FALSE)
huc12_ids <- gsub("huc12_|_raster1985loss.tif", "", huc12_files)
lookup_table <- lookup_table[sapply(lookup_table$huc12s, function(x) any(x %in% huc12_ids)), ]
service_area_indices <- seq_len(nrow(lookup_table))  # Adjust indexing

num_service_areas <- nrow(lookup_table)
cat("Total service areas with wetland loss:", num_service_areas, "\n\n")

process_service_area <- function(i, lookup_table) {
  
  setwd("L:/Wetland Flood Mitigation")
  library(terra)
  
  sa_polygons <- readRDS("Service Areas/ServiceAreas_agg.rds")
  
  output_dir  <- "Wetland Loss/Service Area Mosaics Loss"
  raster_dir  <- "Wetland Loss/Output loss 2021"
  lookup_table_path <- "Service Areas/Service_Area_HUC12_Lookup.rds"
  
  lookup_table <- readRDS(lookup_table_path)
  lookup_table <- lookup_table[sapply(lookup_table$huc12s, function(x) any(x %in% huc12_ids)), ]
  
  sa_id  <- lookup_table$ID[i]
  out_path <- file.path(output_dir, paste0(sa_id, "_lossmosaic.tif"))
  
  if (file.exists(out_path)) {
    message(sprintf("Skipping service area %s (already processed).", sa_id))
    return(NULL)
  }
  
  huc12s <- lookup_table$huc12s[[i]]
  huc12s <- intersect(huc12s, huc12_ids)  # Keep only those with data
  
  if (length(huc12s) == 0) {
    warning(sprintf("Skipping service area %s (no valid HUC12 rasters).", sa_id))
    return(NULL)
  }
  
  message(sprintf("Processing service area %d of %d (ID: %s)", i, num_service_areas, sa_id))
  
  tif_paths <- file.path(raster_dir, paste0("huc12_", huc12s, "_raster1985loss.tif"))
  tif_paths <- tif_paths[file.exists(tif_paths)]  
  
  message(sprintf(" - Reading %d raster files...", length(tif_paths)))
  ras_list <- lapply(tif_paths, rast)
  
  message(" - Mosaicking rasters...")
  mosaic_ras <- do.call(terra::mosaic, ras_list)
  
  message(" - Extracting service area polygon...")
  sa_poly <- sa_polygons[sa_polygons$ID == sa_id, ]
  if (nrow(sa_poly) == 0) {
    warning(sprintf("No matching polygon found for Service Area %s.", sa_id))
    return(NULL)
  }
  
  message(" - Cropping and masking mosaic to service area polygon...")
  mosaic_cropped <- terra::crop(mosaic_ras, sa_poly)
  mosaic_masked <- terra::mask(mosaic_cropped, sa_poly)
  
  message(sprintf(" - Saving mosaic to: %s", out_path))
  terra::writeRaster(mosaic_masked, out_path, overwrite = FALSE)
  
  msg <- sprintf("Mosaic successfully created for service area %s", sa_id)
  message(" - ", msg)
  return(msg)
}

# Retry mechanism
max_retries <- 5  
retry_count <- 0

plan(multisession, workers = 5)

repeat {
  retry_count <- retry_count + 1
  message("Attempt: ", retry_count)
  
  results <- tryCatch({
    future_lapply(
      X = service_area_indices,
      FUN = function(i) process_service_area(i, lookup_table),
      future.seed = TRUE
    )
  }, error = function(e) {
    message("Error encountered: ", e$message)
    return(NULL)
  })
  
  if (!is.null(results) || retry_count >= max_retries) {
    break
  }
}
