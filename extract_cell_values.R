
library(terra)
library(progressr)
library(data.table)

setwd("C:/Users/indumati/Box/Paper2_final")

# --------------------------------------------------------#
#-                      Extraction - SA                  -#
#---------------------------------------------------------#


# Paths
raster_folder <- "Service Area Mosaics"  # Folder where service area TIFFs are stored
output_folder <- "Extractions and Summaries/Extract Service Areas"

# Ensure output directory exists
if (!dir.exists(output_folder)) {
  dir.create(output_folder, recursive = TRUE)
}

# Load service area dataset
sas <- readRDS("Service Areas/ServiceAreas_agg.rds")

# Get list of already processed files
existing_files <- list.files(output_folder, pattern = "^sa_extract_.*\\.rds$", full.names = FALSE)
processed_sa_ids <- gsub("^sa_extract_(.*)\\.rds$", "\\1", existing_files)  # Extract service area IDs

# Identify service areas that still need processing
sas_to_process <- sas[!(sas$ID %in% processed_sa_ids), ]

cat("Total service areas:", nrow(sas), "\n")
cat("Already processed:", length(processed_sa_ids), "\n")
cat("Remaining to process:", nrow(sas_to_process), "\n\n")

# Exit early if everything is processed
if (nrow(sas_to_process) == 0) {
  cat("All service areas have been processed. Exiting...\n")
  quit(save = "no")
}

handlers(global = TRUE)
progressr::with_progress({
  p <- progressr::progressor(along = sas_to_process$ID)
  
  failed_services <- c()  # Store IDs of failed extractions
  
  results <- lapply(1:nrow(sas_to_process), function(i) {
    service_area <- sas_to_process[i, ]
    service_area_name <- service_area$ID
    raster_path <- file.path(raster_folder, paste0(service_area_name, "_mosaic.tif"))
    output_path <- file.path(output_folder, paste0("sa_extract_", service_area_name, ".rds"))
    
    p(sprintf("Processing: %s", service_area_name))
    
    # Check if raster exists
    if (!file.exists(raster_path)) {
      warning(sprintf("Raster not found for %s, skipping...", service_area_name))
      failed_services <<- c(failed_services, service_area_name)
      return(NULL)
    }
    
    # Load raster
    sa_raster <- rast(raster_path)
    
    # Attempt extraction with error handling
    extracted_data <- tryCatch({
      terra::extract(x = sa_raster, y = service_area, cells = TRUE)
    }, error = function(e) {
      warning(sprintf("Extraction failed for %s: %s", service_area_name, e$message))
      failed_services <<- c(failed_services, service_area_name)
      return(NULL)
    })
    
    # Skip if extraction failed
    if (is.null(extracted_data)) return(NULL)
    
    # Attach service area name
    extracted_data$service_area_name <- service_area_name  
    
    # Save the extracted data
    saveRDS(extracted_data, output_path)
    
    return(extracted_data)
  })
  
  # Write failed service area IDs to a log file
  if (length(failed_services) > 0) {
    writeLines(failed_services, error_log)
    message("Some extractions failed. Check ", error_log, " for details.")
  }
})

message("Service area extraction complete. All files saved in: ", output_folder)



# --------------------------------------------------------#
#-                      Extraction - BANK                -#
#---------------------------------------------------------#

# Paths
raster_folder <- "Service Area Mosaics"  # Folder where service area TIFFs are stored
output_folder <- "Extractions and Summaries/Extract Banks"

# Load datasets
banks <- readRDS("Bank Footprints/footprints_and_buffers.rds")
sas <- readRDS("Service Areas/ServiceAreas_agg.rds")
  # rows_to_update <- which(sas$ID %in% c("Three_Lakes_Regional_MB_FDOT_", "FP_L_Everglades_Phase_II_MB"))
  # sas$ID[rows_to_update] <- c("FP_AndL_Everglades_Phase_II_MB", "Three_Lakes_Regional_MB_FDOT")


# Function to sanitize names
sanitize_name <- function(name) {
  name <- gsub("[^A-Za-z0-9_-]", "_", name)  # Replace non-alphanumeric characters with "_"
  name <- gsub("_+", "_", name)  # Replace multiple "_" with a single "_"
  name <- gsub("^_|_$", "", name)  # Remove leading or trailing "_"
  return(name)
}

# Apply name sanitization
banks$Name <- sapply(banks$Name, sanitize_name)
#saveRDS(banks, "Footprints/footprints_and_buffers.rds")

# Filter banks to match those in service areas
banks <- banks[banks$Name %in% sas$ID, ]

# Create output directory if it doesn't exist
if (!dir.exists(output_folder)) {
  dir.create(output_folder, recursive = TRUE)
}

# Progress bar setup
handlers(global = TRUE)
progressr::with_progress({
  p <- progressr::progressor(along = banks$Name)
  
  results <- lapply(1:nrow(banks), function(i) {
    test_bank <- banks[i, ]
    test_bank_name <- test_bank$Name
    raster_path <- file.path(raster_folder, paste0(test_bank_name, "_mosaic.tif"))
    
    p(sprintf("Processing: %s", test_bank_name))
    
    # Check if raster exists
    if (!file.exists(raster_path)) {
      warning(sprintf("Raster not found for %s, skipping...", test_bank_name))
      return(NULL)
    }
    
    # Load raster
    sa_raster <- rast(raster_path)
    
    # Extract values
    extracted_data <- terra::extract(
      x = sa_raster,
      y = test_bank,  # Extract for the single bank polygon
      cells = TRUE
    )
    
    # Attach bank name
    extracted_data$bank_name <- test_bank_name  
    
    # Save the extracted data
    saveRDS(extracted_data, file.path(output_folder, paste0("bank_extract_", test_bank_name, ".rds")))
    
    return(extracted_data)
  })
})

message("Bank extraction complete. All files saved in: ", output_folder)


# --------------------------------------------------------#
#-            Extractions for Outlier Banks              -#
#---------------------------------------------------------#

# Paths
raster_folder <- "CONUS\\Service Area Mosaics 2021"  # Folder where service area TIFFs are stored
output_folder <- "Extractions and Summaries/Extract Banks"

# Load datasets
outlier_banks <- readRDS("Bank Footprints/outlier_banks.rds")

sas <- readRDS("Service Areas/ServiceAreas_agg.rds")

# Function to sanitize names
sanitize_name <- function(name) {
  name <- gsub("[^A-Za-z0-9_-]", "_", name)  # Replace non-alphanumeric characters with "_"
  name <- gsub("_+", "_", name)  # Replace multiple "_" with a single "_"
  name <- gsub("^_|_$", "", name)  # Remove leading or trailing "_"
  return(name)
}

# Apply name sanitization
outlier_banks$Name <- sapply(outlier_banks$Name, sanitize_name)
#saveRDS(banks, "Footprints/footprints_and_buffers.rds")

# Filter banks to match those in service areas
outlier_banks <- outlier_banks[outlier_banks$Name %in% sas$ID, ]

# Create output directory if it doesn't exist
if (!dir.exists(output_folder)) {
  dir.create(output_folder, recursive = TRUE)
}

# Progress bar setup
handlers(global = TRUE)
progressr::with_progress({
  p <- progressr::progressor(along = outlier_banks$Name)
  
  results <- lapply(1:nrow(outlier_banks), function(i) {
    test_bank <- outlier_banks[i, ]
    test_bank_name <- test_bank$Name
    raster_path <- file.path(raster_folder, paste0(test_bank_name, "_mosaic.tif"))
    
    p(sprintf("Processing: %s", test_bank_name))
    
    # Check if raster exists
    if (!file.exists(raster_path)) {
      warning(sprintf("Raster not found for %s, skipping...", test_bank_name))
      return(NULL)
    }
    
    # Load raster
    sa_raster <- rast(raster_path)
    
    # Extract values
    extracted_data <- terra::extract(
      x = sa_raster,
      y = test_bank,  # Extract for the single bank polygon
      cells = TRUE
    )
    
    # Attach bank name
    extracted_data$bank_name <- test_bank_name  
    
    # Save the extracted data
    saveRDS(extracted_data, file.path(output_folder, paste0("bank_extract_", test_bank_name, ".rds")))
    
    return(extracted_data)
  })
})

message("Bank extraction complete. All outlier bank files saved in: ", output_folder)