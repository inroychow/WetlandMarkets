# --------------------------------------------------------
#                   Packages & Setup
# --------------------------------------------------------

library(terra)
library(sf)
library(dplyr)
library(future.apply)
library(fst)
library(progressr)
library(data.table)

# Base directory 
datadir <- Sys.getenv("DATA_DIR", "C:/Users/indumati/Box/Paper2_final")

# Helper to build absolute paths quickly
pth <- function(...) file.path(datadir, ...)

# --------------------------------------------------------
#             Extraction: Service Areas (SA)
# --------------------------------------------------------

raster_folder <- pth("Service Area Mosaics")
output_folder <- pth("Extractions and Summaries", "Extract Service Areas")

if (!dir.exists(output_folder)) dir.create(output_folder, recursive = TRUE)

sas <- readRDS(pth("Service Areas", "ServiceAreas_agg.rds"))

existing_files <- list.files(output_folder, pattern = "^sa_extract_.*\\.rds$", full.names = FALSE)
processed_sa_ids <- gsub("^sa_extract_(.*)\\.rds$", "\\1", existing_files)
sas_to_process <- sas[!(sas$ID %in% processed_sa_ids), ]

cat("Total service areas:", nrow(sas), "\n")
cat("Already processed:", length(processed_sa_ids), "\n")
cat("Remaining to process:", nrow(sas_to_process), "\n\n")

if (nrow(sas_to_process) == 0) {
  cat("All service areas have been processed. Exiting...\n")
  quit(save = "no")
}

handlers(global = TRUE)
progressr::with_progress({
  p <- progressr::progressor(along = sas_to_process$ID)
  failed_services <- c()
  
  results <- lapply(1:nrow(sas_to_process), function(i) {
    service_area <- sas_to_process[i, ]
    sa_id <- service_area$ID
    raster_path <- file.path(raster_folder, paste0(sa_id, "_mosaic.tif"))
    output_path <- file.path(output_folder, paste0("sa_extract_", sa_id, ".rds"))
    
    p(sprintf("Processing: %s", sa_id))
    
    if (!file.exists(raster_path)) {
      warning(sprintf("Raster not found for %s, skipping...", sa_id))
      failed_services <<- c(failed_services, sa_id)
      return(NULL)
    }
    
    sa_raster <- rast(raster_path)
    
    extracted_data <- tryCatch({
      terra::extract(sa_raster, service_area, cells = TRUE)
    }, error = function(e) {
      warning(sprintf("Extraction failed for %s: %s", sa_id, e$message))
      failed_services <<- c(failed_services, sa_id)
      return(NULL)
    })
    
    if (is.null(extracted_data)) return(NULL)
    
    extracted_data$service_area_name <- sa_id
    saveRDS(extracted_data, output_path)
    
    return(extracted_data)
  })
  
  if (length(failed_services) > 0) {
    writeLines(failed_services, pth("Logs", "extraction_failures_sa.txt"))
    message("Some extractions failed. Check extraction_failures_sa.txt for details.")
  }
})

message("Service area extraction complete. All files saved in: ", output_folder)

# --------------------------------------------------------
#             Extraction: Banks (Standard)
# --------------------------------------------------------

raster_folder <- pth("Service Area Mosaics")
output_folder <- pth("Extractions and Summaries", "Extract Banks")

banks <- readRDS(pth("Bank Footprints", "footprints_and_buffers.rds"))
sas <- readRDS(pth("Service Areas", "ServiceAreas_agg.rds"))

sanitize_name <- function(name) {
  name <- gsub("[^A-Za-z0-9_-]", "_", name)
  name <- gsub("_+", "_", name)
  name <- gsub("^_|_$", "", name)
  return(name)
}

banks$Name <- sapply(banks$Name, sanitize_name)
banks <- banks[banks$Name %in% sas$ID, ]

if (!dir.exists(output_folder)) dir.create(output_folder, recursive = TRUE)

handlers(global = TRUE)
progressr::with_progress({
  p <- progressr::progressor(along = banks$Name)
  
  results <- lapply(1:nrow(banks), function(i) {
    bank <- banks[i, ]
    bank_name <- bank$Name
    raster_path <- file.path(raster_folder, paste0(bank_name, "_mosaic.tif"))
    
    p(sprintf("Processing: %s", bank_name))
    
    if (!file.exists(raster_path)) {
      warning(sprintf("Raster not found for %s, skipping...", bank_name))
      return(NULL)
    }
    
    sa_raster <- rast(raster_path)
    extracted_data <- terra::extract(sa_raster, bank, cells = TRUE)
    extracted_data$bank_name <- bank_name
    
    saveRDS(extracted_data, file.path(output_folder, paste0("bank_extract_", bank_name, ".rds")))
    return(extracted_data)
  })
})

message("Bank extraction complete. All files saved in: ", output_folder)

# --------------------------------------------------------
#             Extraction: Outlier Banks
# --------------------------------------------------------

raster_folder <- pth("CONUS", "Service Area Mosaics 2021")
output_folder <- pth("Extractions and Summaries", "Extract Banks")

outlier_banks <- readRDS(pth("Bank Footprints", "outlier_banks.rds"))
sas <- readRDS(pth("Service Areas", "ServiceAreas_agg.rds"))

outlier_banks$Name <- sapply(outlier_banks$Name, sanitize_name)
outlier_banks <- outlier_banks[outlier_banks$Name %in% sas$ID, ]

if (!dir.exists(output_folder)) dir.create(output_folder, recursive = TRUE)

handlers(global = TRUE)
progressr::with_progress({
  p <- progressr::progressor(along = outlier_banks$Name)
  
  results <- lapply(1:nrow(outlier_banks), function(i) {
    bank <- outlier_banks[i, ]
    bank_name <- bank$Name
    raster_path <- file.path(raster_folder, paste0(bank_name, "_mosaic.tif"))
    
    p(sprintf("Processing: %s", bank_name))
    
    if (!file.exists(raster_path)) {
      warning(sprintf("Raster not found for %s, skipping...", bank_name))
      return(NULL)
    }
    
    sa_raster <- rast(raster_path)
    extracted_data <- terra::extract(sa_raster, bank, cells = TRUE)
    extracted_data$bank_name <- bank_name
    
    saveRDS(extracted_data, file.path(output_folder, paste0("bank_extract_", bank_name, ".rds")))
    return(extracted_data)
  })
})

message("Outlier bank extraction complete. All files saved in: ", output_folder)

# --------------------------------------------------------
#         Combine All Bank Extraction Files
# --------------------------------------------------------

bank_extract_folder <- pth("Extractions and Summaries", "Extract Banks")

bank_files <- list.files(bank_extract_folder, pattern = "^bank_extract_.*\\.rds$", full.names = TRUE)

bank_extract_combined <- bind_rows(lapply(bank_files, readRDS))

saveRDS(bank_extract_combined, file = file.path(bank_extract_folder, "bank_extractions.rds"))

message("All individual bank extractions combined into 'bank_extractions.rds'.")

# --------------------------------------------------------
#      Extraction: Wetland Loss - Service Areas (SA)
# --------------------------------------------------------

# Base folders for wetland loss raster mosaics
raster_folder <- pth("Service Area Mosaics", "Loss Mosaics")
output_folder_loss <- pth("Extractions and Summaries", "Extract Loss")

dir.create(output_folder_loss, recursive = TRUE, showWarnings = FALSE)

sas <- readRDS(pth("Service Areas", "ServiceAreas_agg.rds"))

existing_files <- list.files(output_folder_loss, pattern = "^loss_extract_.*\\.rds$", full.names = FALSE)
processed_sa_ids <- gsub("^loss_extract_(.*)\\.rds$", "\\1", existing_files)
sas_to_process <- sas[!(sas$ID %in% processed_sa_ids), ]

cat("Wetland Loss: Total service areas:", nrow(sas), "\n")
cat("Already processed:", length(processed_sa_ids), "\n")
cat("Remaining to process:", nrow(sas_to_process), "\n\n")

handlers(global = TRUE)
progressr::with_progress({
  p <- progressor(along = sas_to_process$ID)
  failed_sas <- c()
  
  lapply(1:nrow(sas_to_process), function(i) {
    sa <- sas_to_process[i, ]
    sa_id <- sa$ID
    raster_path <- file.path(raster_folder, paste0(sa_id, "_loss1985mosaic.tif"))
    output_path <- file.path(output_folder_loss, paste0("loss_extract_", sa_id, ".rds"))
    
    p(sprintf("Processing SA (loss): %s", sa_id))
    
    if (!file.exists(raster_path)) {
      warning(sprintf("Raster not found for SA %s", sa_id))
      failed_sas <<- c(failed_sas, sa_id)
      return(NULL)
    }
    
    sa_raster <- rast(raster_path)
    extracted_data <- tryCatch({
      terra::extract(sa_raster, sa, cells = TRUE)
    }, error = function(e) {
      warning(sprintf("Extraction failed for SA %s: %s", sa_id, e$message))
      failed_sas <<- c(failed_sas, sa_id)
      return(NULL)
    })
    
    if (is.null(extracted_data)) return(NULL)
    extracted_data$service_area_name <- sa_id
    saveRDS(extracted_data, output_path)
    return(extracted_data)
  })
  
  if (length(failed_sas) > 0) {
    writeLines(failed_sas, pth("Logs", "extraction_failures_sa_loss.txt"))
    message("Some loss extractions failed. Check 'extraction_failures_sa_loss.txt'.")
  }
})

message("Wetland loss extractioncomplete.")