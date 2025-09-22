#--------------------------------------------------------------
#  Mosaic HUC12 rasters by Service Area – for 2021 and loss
#--------------------------------------------------------------
# 
#         Set DATA_DIR env‑var (or edit datadir) and run.
#         Call `mosaic_service_areas()` once for each raster set (loss or 2021)
#         you need to mosaic 
#--------------------------------------------------------------

# Packages ----------------------------------------------------
library(terra)
library(sf)
library(dplyr)
library(future.apply)
library(fst)

# Base directory (override with env‑var) ----------------------
datadir <- Sys.getenv("DATA_DIR", "C:/Users/indumati/Box/Paper2_final")

# Helper to build absolute paths quickly
pth <- function(...) file.path(datadir, ...)

#--------------------------------------------------------------
#  Generic mosaicking function
#--------------------------------------------------------------

mosaic_service_areas <- function(lookup_table_path,
                                 sa_polygons_path,
                                 raster_dir,
                                 raster_pattern,          # regex with one capture group for HUC12 id
                                 raster_filename_tmpl,    # sprintf‑style template with %s for HUC12 id
                                 output_dir,
                                 mosaic_suffix = "_mosaic.tif",
                                 workers      = 5,
                                 skip_existing = TRUE) {
  # ---- Input validation ------------------------------------
  stopifnot(dir.exists(dirname(lookup_table_path)),
            dir.exists(dirname(sa_polygons_path)),
            dir.exists(raster_dir))
  if (!dir.exists(output_dir)) dir.create(output_dir, recursive = TRUE)
  
  # ---- Read static objects ---------------------------------
  lookup_table <- readRDS(lookup_table_path)
  sa_polygons  <- readRDS(sa_polygons_path)
  
  # ---- Detect available rasters ----------------------------
  raster_files <- list.files(raster_dir, pattern = "tif$", full.names = FALSE)
  huc12_ids    <- sub(raster_pattern, "\\1", raster_files, perl = TRUE)
  available    <- huc12_ids != raster_files      # keep only those that matched
  huc12_ids    <- huc12_ids[available]
  raster_files <- raster_files[available]
  
  # ---- Filter lookup to service areas that actually have data
  lookup_table <- lookup_table[sapply(lookup_table$huc12s,
                                      function(x) any(x %in% huc12_ids)), ]
  n_sa <- nrow(lookup_table)
  message("Total service areas to process: ", n_sa)
  if (n_sa == 0) {
    warning("No service areas with matching rasters found – nothing to do.")
    return(invisible(NULL))
  }
  
  # ---- Internal worker function ----------------------------
  process_sa <- function(idx) {
    sa_id   <- lookup_table$ID[idx]
    outfile <- file.path(output_dir, paste0(sa_id, mosaic_suffix))
    
    if (skip_existing && file.exists(outfile)) {
      message("[", sa_id, "] skipped (exists)")
      return(sa_id)
    }
    
    huc12s <- intersect(lookup_table$huc12s[[idx]], huc12_ids)
    if (length(huc12s) == 0) {
      warning("[", sa_id, "] skipped – no rasters available")
      return(NA)
    }
    
    # Build full paths for this SA
    tif_paths <- file.path(raster_dir, sprintf(raster_filename_tmpl, huc12s))
    tif_paths <- tif_paths[file.exists(tif_paths)]
    
    # Read & mosaic
    ras_list  <- lapply(tif_paths, rast)
    mosaic    <- do.call(terra::mosaic, ras_list)
    
    # Crop & mask to SA polygon
    sa_poly   <- sa_polygons[sa_polygons$ID == sa_id, ]
    if (nrow(sa_poly) == 0) {
      warning("[", sa_id, "] no matching polygon – skipped")
      return(NA)
    }
    mosaic_sa <- terra::mask(terra::crop(mosaic, sa_poly), sa_poly)
    
    # Save
    terra::writeRaster(mosaic_sa, outfile, overwrite = TRUE)
    message("[", sa_id, "] mosaic written → ", basename(outfile))
    sa_id
  }
  
  # ---- Parallel execution ----------------------------------
  plan(multisession, workers = workers)
  on.exit(plan(sequential), add = TRUE)
  results <- future_lapply(seq_len(n_sa), process_sa, future.seed = TRUE)
  invisible(results)
}

#--------------------------------------------------------------
#  EXAMPLE CALLS                        -----------------------
#--------------------------------------------------------------

# 1) Classified HUC12 rasters --------------------------------
# mosaic_service_areas(
#   lookup_table_path   = pth("Service Areas", "Service_Area_HUC12_Lookup.rds"),
#   sa_polygons_path    = pth("Service Areas", "ServiceAreas_agg.rds"),
#   raster_dir          = pth("Classified HUC12 Rasters"),
#   raster_pattern      = "huc12_(.*)_raster.tif$",
#   raster_filename_tmpl= "huc12_%s_raster.tif",
#   output_dir          = pth("Mosaics"),
#   workers             = 5)

# 2) Wetland loss (1985‑2021) --------------------------------
# mosaic_service_areas(
#   lookup_table_path   = pth("Service Areas", "Service_Area_HUC12_Lookup.rds"),
#   sa_polygons_path    = pth("Service Areas", "ServiceAreas_agg.rds"),
#   raster_dir          = pth("Wetland Loss", "Classified HUC12 Rasters - Loss"),
#   raster_pattern      = "huc12_(.*)_raster1985loss.tif$",
#   raster_filename_tmpl= "huc12_%s_raster1985loss.tif",
#   output_dir          = pth("Service Area Mosaics", "Loss Mosaics"),
#   workers             = 5)

