library(terra)
library(sf)
library(dplyr)
library(future.apply)
library(fst) 

setwd("C:/Users/indumati/Box/Paper2_final")


#-----------------------------------------------#
#           To mosaic sequentially ...          #
#-----------------------------------------------#

output_dir <- "Mosaics"
raster_dir <- "Classified HUC12 Rasters"
lookup_file <- "Service Areas/Service_Area_HUC12_Lookup.rds"

sa_polygons <- readRDS("Service Areas/ServiceAreas_agg.rds")  
lookup_table <- readRDS("Service Areas/Service_Area_HUC12_Lookup.rds")

# Load necessary library
library(terra)

# Filter out already processed service areas
already_processed <- list.files(output_dir, pattern = "_mosaic.tif$", full.names = FALSE)
processed_ids <- gsub("_mosaic.tif", "", already_processed)
lookup_table <- lookup_table[!(lookup_table$ID %in% processed_ids), ]

# Print the number of rows in lookup_table
print(paste("Total service areas to process:", nrow(lookup_table)))

# Run sequentially with lapply
results <- lapply(seq_len(nrow(lookup_table)), function(i) {
  # Print progress
  print(paste("Processing service area", i, "of", nrow(lookup_table)))
  
  # Extract service area ID and HUC12 codes
  sa_id <- lookup_table$ID[i]
  print(paste("Service Area ID:", sa_id))
  
  huc12s <- lookup_table$huc12s[[i]]
  print(paste("Number of HUC12s:", length(huc12s)))
  
  # Build full paths to the TIF files for these HUC12s
  tif_paths <- file.path(raster_dir, paste0("huc12_", huc12s, "_raster.tif"))
  tif_paths <- tif_paths[file.exists(tif_paths)]  # keep only existing files
  print(paste("Valid raster files found:", length(tif_paths)))
  
  # Skip if no rasters are found
  if (length(tif_paths) == 0) {
    print(paste("Skipping service area", sa_id, "due to no valid rasters."))
    return(NULL)
  }
  
  # Read rasters and create mosaic
  print(paste("Reading", length(tif_paths), "raster files..."))
  ras_list <- lapply(tif_paths, function(path) {
    print(paste("Reading raster:", path))
    rast(path)
  })
  
  print("Mosaicking rasters...")
  mosaic_ras <- do.call(terra::mosaic, ras_list)
  print("Mosaic complete.")
  
  # Get the polygon for this service area
  print("Extracting service area polygon...")
  sa_poly <- sa_polygons[sa_polygons$ID == sa_id, ]
  
  # Ensure sa_poly is not empty
  if (nrow(sa_poly) == 0) {
    print(paste("Warning: No matching polygon found for Service Area", sa_id))
    return(NULL)
  }
  
  print("Cropping mosaic to service area polygon...")
  mosaic_cropped <- crop(mosaic_ras, sa_poly)
  
  print("Masking mosaic to service area polygon...")
  mosaic_masked <- mask(mosaic_cropped, sa_poly)
  
  # Define output path
  out_path <- file.path(output_dir, paste0(sa_id, "_mosaic.tif"))
  print(paste("Saving mosaic to:", out_path))
  
  # Write out the mosaic
  writeRaster(mosaic_masked, out_path, overwrite = TRUE)
  
  print(paste("Mosaic successfully created for service area", sa_id))
  
  # Return a message to track progress
  return(paste("Mosaic created for service area", sa_id))
})

print("Processing complete!")


#-----------------------------------------------#
#           To mosaic in parallel ...           #
#-----------------------------------------------#


# Paths
output_dir        <- "Mosaics"
raster_dir        <- "Classified HUC12 Rasters"

# Load SA polygons 
sa_polygons       <- readRDS("Service Areas/ServiceAreas_agg.rds")

# Load the lookup table and subset it correctly to the SAs you want to run
lookup_table <- readRDS("Service Areas/Service_Area_HUC12_Lookup.rds")
lookup_table <- lookup_table[500:1863, ]  # Subset before use
service_area_indices <- seq(500, 1863)  # Explicit indices

num_service_areas <- nrow(lookup_table)
cat("Total service areas:", num_service_areas, "\n\n")

#Load ALL objects in function for parallel processing.
process_service_area <- function(i, lookup_table) {
  
  # Load libraries inside the function for parallel execution
  setwd("L:/Wetland Flood Mitigation")
  library(terra)
  sa_polygons       <- readRDS("Service Areas/ServiceAreas_agg.rds")
  
  output_dir  <- "Mosaics"
  raster_dir  <- "Classified HUC12 Rasters"
  lookup_table_path <- "Service Areas/Service_Area_HUC12_Lookup.rds"
  
  lookup_table <- readRDS(lookup_table_path)
  lookup_table <- lookup_table[500:1863, ]  # Subset before use
  
  # Extract service area ID and HUC12 codes
  sa_id  <- lookup_table$ID[i]
  
  # Define output path
  out_path <- file.path(output_dir, paste0(sa_id, "_mosaic.tif"))
  
  # Check if the file already exists
  if (file.exists(out_path)) {
    message(sprintf("Skipping service area %s (already processed).", sa_id))
    return(NULL)  # Skip processing if already done
  }
  
  huc12s <- lookup_table$huc12s[[i]]
  
  message(sprintf("Processing service area %d of %d (ID: %s)", i, num_service_areas, sa_id))
  
  # Build full paths to the TIF files for these HUC12s
  tif_paths <- file.path(raster_dir, paste0("huc12_", huc12s, "_raster.tif"))
  tif_paths <- tif_paths[file.exists(tif_paths)]  # Keep only existing files
  
  if (length(tif_paths) == 0) {
    warning(sprintf("Skipping service area %s (no valid rasters).", sa_id))
    return(NULL)
  }
  
  # Read rasters and mosaic
  message(sprintf(" - Reading %d raster files...", length(tif_paths)))
  ras_list <- lapply(tif_paths, rast)
  
  message(" - Mosaicking rasters...")
  mosaic_ras <- do.call(terra::mosaic, ras_list)
  
  message(" - Extract the polygon for this service area...")
  
  # Extract the polygon for this service area
  sa_poly <- sa_polygons[sa_polygons$ID == sa_id, ]
  if (nrow(sa_poly) == 0) {
    warning(sprintf("No matching polygon found for Service Area %s.", sa_id))
    return(NULL)
  }
  
  message(" - Cropping and masking mosaic to service area polygon...")
  mosaic_cropped <- terra::crop(mosaic_ras, sa_poly)
  mosaic_masked <- terra::mask(mosaic_cropped, sa_poly)
  
  message(sprintf(" - Saving mosaic to: %s", out_path))
  
  # Write mosaic to disk
  terra::writeRaster(mosaic_masked, out_path, overwrite = FALSE)
  
  msg <- sprintf("Mosaic successfully created for service area %s", sa_id)
  message(" - ", msg)
  return(msg)
}

# Set maximum retries in case of failure
max_retries <- 5  
retry_count <- 0

plan(multisession, workers = 10)

repeat {
  retry_count <- retry_count + 1
  message("Attempt: ", retry_count)
  
  results <- tryCatch({
    future_lapply(
      X = service_area_indices,  # Use correct indexing
      FUN = function(i) process_service_area(i, lookup_table),  # Pass lookup_table
      future.seed = TRUE
    )
  }, error = function(e) {
    message("Error encountered: ", e$message)
    return(NULL)  # Return NULL on error so we can check if it failed
  })
  
  if (!is.null(results) || retry_count >= max_retries) {
    break  # Exit loop if successful or max retries reached
  }
}

plan(sequential)