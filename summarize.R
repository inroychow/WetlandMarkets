library(terra)
library(progressr)
library(sf)
library(dplyr)
library(future.apply)

setwd("C:/Users/indumati/Box/Paper2_final")

# ---------------------------------------------- #
# ------------- Summarize banks -----------------# 
# ---------------------------------------------- #

example_rast = rast("Mosaics/MT-Running_Colter_Ranch_mosaic.tif")

banks = vect("Bank Footprints/footprints_and_buffers.rds")
banks = project(banks, crs(final_multilayer))

# ---------------------------------------------------------------------------------------

# Paths to spatial datasets
huc12_shp <- "Variables/huc12_boundaries.gpkg"
sa_shp <- "Service Areas/ServiceAreas_agg.rds"
output_folder <- "Service Areas/Temp lookup"


## Summarize bank -----------------------------------------------------------------------
bank_extract = readRDS("all_bank_extractions.rds")
# colnames(bank_extract)[2:4] <- c("housing units", "housing value", "population")
bank_summary <- bank_extract %>%
  # Exclude rows with no raster data
  group_by(bank_name) %>%
  summarise(
    sum_housing_units = sum(`housing units`, na.rm = TRUE),
    mean_housing_units = mean(`housing units`, na.rm = TRUE),
    median_housing_units = median(`housing units`, na.rm = TRUE),
    min_housing_units = min(`housing units`, na.rm = TRUE),
    max_housing_units = max(`housing units`, na.rm = TRUE),
    sd_housing_units = sd(`housing units`, na.rm = TRUE),
    q25_housing_units = quantile(`housing units`, 0.25, na.rm = TRUE),
    q75_housing_units = quantile(`housing units`, 0.75, na.rm = TRUE),
    
    sum_housing_value = sum(`housing value`, na.rm = TRUE),
    mean_housing_value = mean(`housing value`, na.rm = TRUE),
    median_housing_value = median(`housing value`, na.rm = TRUE),
    min_housing_value = min(`housing value`, na.rm = TRUE),
    max_housing_value = max(`housing value`, na.rm = TRUE),
    sd_housing_value = sd(`housing value`, na.rm = TRUE),
    q25_housing_value = quantile(`housing value`, 0.25, na.rm = TRUE),
    q75_housing_value = quantile(`housing value`, 0.75, na.rm = TRUE),
    
    sum_population = sum(population, na.rm = TRUE),
    mean_population = mean(population, na.rm = TRUE),
    median_population = median(population, na.rm = TRUE),
    min_population = min(population, na.rm = TRUE),
    max_population = max(population, na.rm = TRUE),
    sd_population = sd(population, na.rm = TRUE),
    q25_population = quantile(population, 0.25, na.rm = TRUE),
    q75_population = quantile(population, 0.75, na.rm = TRUE),
    
    cell_count = n(),  # Total number of cells
    .groups = "drop"
  )

saveRDS(bank_summary, "Extractions and Summaries/Summaries Banks/bank_summary.rds")


# ---------------------------------------------- #
# ------------- Summarize SAs - -----------------# 
# ---------------------------------------------- #

SAs = readRDS("Service Areas/ServiceAreas_agg.rds")


#------------------------------------------------------------------------------
# 1. Setup paths
#------------------------------------------------------------------------------
extracted_dir <- "Extractions/Extract Service Areas"             # Where the .rds extracted data live
output_dir_sa <- "Extractions/Summaries Service Areas"  # Where to save the summary files
dir.create(output_dir_sa, recursive = TRUE, showWarnings = FALSE)

# List all RDS extraction files with a pattern "sa_extract_<SA>.rds"
rds_files <- list.files(
  extracted_dir, 
  pattern = "^sa_extract_.*\\.rds$", 
  full.names = TRUE
)

# Identify already processed service areas
existing_summaries <- list.files(
  output_dir_sa, 
  pattern = "^summary_.*\\.rds$", 
  full.names = FALSE
)
processed_sa_names <- gsub("^summary_(.*)\\.rds$", "\\1", existing_summaries)  # Extract names without "summary_" prefix

#------------------------------------------------------------------------------
# 2. Define helper functions to do the summarizing and ECDF
#------------------------------------------------------------------------------
summarize_extraction <- function(df, sa_name) {
  # Read bank extract
  bank_extract <- readRDS("Extractions/Extract Banks/all_bank_extractions.rds")
  
  # Subset to only the banks associated with this service area
  bank_extract_sa <- bank_extract[bank_extract$bank_name == sa_name, ]
  
  # Extract bank cell IDs
  bank_cell_ids <- bank_extract_sa$cell
  
  # Remove cells that appear in bank_extract using indexing
  df <- df[!df$cell %in% bank_cell_ids, ]
  
  sa_summary <- df %>%
    summarise(
      sum_housing_units    = sum(housing_units, na.rm = TRUE),
      mean_housing_units   = mean(housing_units, na.rm = TRUE),
      median_housing_units = median(housing_units, na.rm = TRUE),
      min_housing_units    = min(housing_units, na.rm = TRUE),
      max_housing_units    = max(housing_units, na.rm = TRUE),
      sd_housing_units     = sd(housing_units, na.rm = TRUE),
      q25_housing_units    = quantile(housing_units, 0.25, na.rm = TRUE),
      q75_housing_units    = quantile(housing_units, 0.75, na.rm = TRUE),
      
      sum_housing_value    = sum(housing_value, na.rm = TRUE),
      mean_housing_value   = mean(housing_value, na.rm = TRUE),
      median_housing_value = median(housing_value, na.rm = TRUE),
      min_housing_value    = min(housing_value, na.rm = TRUE),
      max_housing_value    = max(housing_value, na.rm = TRUE),
      sd_housing_value     = sd(housing_value, na.rm = TRUE),
      q25_housing_value    = quantile(housing_value, 0.25, na.rm = TRUE),
      q75_housing_value    = quantile(housing_value, 0.75, na.rm = TRUE),
      
      sum_population       = sum(population, na.rm = TRUE),
      mean_population      = mean(population, na.rm = TRUE),
      median_population    = median(population, na.rm = TRUE),
      min_population       = min(population, na.rm = TRUE),
      max_population       = max(population, na.rm = TRUE),
      sd_population        = sd(population, na.rm = TRUE),
      q25_population       = quantile(population, 0.25, na.rm = TRUE),
      q75_population       = quantile(population, 0.75, na.rm = TRUE),
      
      cell_count           = n()
    ) %>%
    mutate(sa_name = sa_name)
  
  return(sa_summary)
}

compute_ecdf_data <- function(df, sa_name) {
  
  # Read bank extract
  bank_extract <- readRDS("Extractions/Extract Banks 2021/all_bank_extractions.rds")
  
  # Subset to only the banks associated with this service area
  bank_extract_sa <- bank_extract[bank_extract$bank_name == sa_name, ]
  
  # Extract bank cell IDs
  bank_cell_ids <- bank_extract_sa$cell
  
  # Remove cells that appear in bank_extract using indexing
  df <- df[!df$cell %in% bank_cell_ids, ]
  
  vars <- c("housing_units", "housing_value", "population")
  
  ecdf_list <- list()
  
  for (var in vars) {
    if (!var %in% colnames(df)) next
    valid_values <- na.omit(df[[var]])
    if (length(valid_values) == 0) {
      ecdf_list[[var]] <- data.frame(value = numeric(0), ecdf = numeric(0))
    } else {
      f_ecdf <- ecdf(valid_values)
      vals   <- sort(unique(valid_values))
      ecdf_list[[var]] <- data.frame(
        value = vals,
        ecdf  = f_ecdf(vals)
      )
    }
  }
  return(ecdf_list)
}

#------------------------------------------------------------------------------
# 3. Loop over each .rds file, summarize, and save
#------------------------------------------------------------------------------
for (file_path in rds_files) {
  
  # Extract the service area name from the filename
  base_name <- basename(file_path) 
  sa_name   <- gsub("^sa_extract_(.*)\\.rds$", "\\1", base_name)
  
  # **Skip if already processed**
  if (sa_name %in% processed_sa_names) {
    cat("\nSkipping already processed:", sa_name, "\n")
    next
  }
  
  cat("\n---------------------------------------------------\n")
  cat("Processing Service Area:", sa_name, "\n")
  
  # 3a) Try reading the extracted data with error handling
  df <- tryCatch({
    readRDS(file_path)
  }, error = function(e) {
    cat("  ERROR reading file:", file_path, "\n  Message:", e$message, "\n  Skipping.\n")
    return(NULL)  # Return NULL to indicate failure
  })
  
  # Skip this service area if reading failed
  if (is.null(df)) next
  
  # 3b) Rename columns (if needed)
  if (ncol(df) >= 4) {
    colnames(df)[2:4] <- c("housing_units", "housing_value", "population")
  }
  
  # 3c) Summarize
  sa_summary <- summarize_extraction(df, sa_name)
  
  # 3d) Compute ECDF data
  ecdf_data  <- compute_ecdf_data(df, sa_name)
  
  # 3e) Save summary as RDS + CSV
  summary_rds_path <- file.path(output_dir_sa, paste0("summary_", sa_name, ".rds"))
  saveRDS(sa_summary, summary_rds_path)
  
  
  
  cat("  Summary saved:", summary_rds_path, "\n")
  
  # 3f) Save ECDF data to RDS
  ecdf_rds_path <- file.path(output_dir_sa, paste0("ecdf_", sa_name, ".rds"))
  saveRDS(ecdf_data, ecdf_rds_path)
  
  cat("  ECDF data saved:", ecdf_rds_path, "\n")
}

##-----------------------------------------##
##          COMBINE INTO ONE SA DF         ##
##-----------------------------------------##

# Define output directory
output_dir <- "Extractions and Summaries/Summaries Service Areas/"

# List all RDS files
rds_files <- list.files(output_dir, pattern = "^summary_.*\\.rds$", full.names = TRUE)

# Read each RDS file into a list
extracted_list <- lapply(rds_files, function(file) {
  df <- readRDS(file)
  
  return(df)
})


# Combine all into a single dataframe
combined_df <- bind_rows(extracted_list, .id = "source")

# Optional: Extract bank name from filename for tracking
combined_df$sa_id <- gsub("summary_(.*)\\.rds", "\\1", basename(rds_files)[as.integer(combined_df$source)])

# Remove temporary 'source' column
combined_df <- select(combined_df, -source)
combined_df <- select(combined_df, -sa_name)

print(head(combined_df))


# Change a faulty bank name
sa_summary$sa_id[sa_summary$sa_id == "Three_Lakes_Regional_MB_FDOT_"] <- "Three_Lakes_Regional_MB_FDOT"
# loss_summary = readRDS("Extractions and Summaries/Summaries Loss/loss_summary.rds")
# loss_summary$sa_id[loss_summary$sa_id == "Three_Lakes_Regional_MB_FDOT_"] <- "Three_Lakes_Regional_MB_FDOT"

# Save the combined dataframe
saveRDS(combined_df, file.path(output_dir, "sa_summary.rds"))


# -------------------

# Combine into one df

bank_summary = readRDS("Extractions and Summaries/Summaries Banks/bank_summary.rds")
sa_summary = readRDS("Extractions and Summaries/Summaries Service Areas/sa_summary.rds")
loss_summary = readRDS("Extractions and Summaries/Summaries Loss/loss_summary.rds")

bank_summary = bank_summary %>% 
  dplyr::rename(id = bank_name,
                sum_hu_bank = sum_housing_units,
                mean_hu_bank = mean_housing_units,
                median_hu_bank = median_housing_units,
                max_hu_bank = max_housing_units) %>% 
  dplyr::select(id, 
                sum_hu_bank,
         mean_hu_bank,
         median_hu_bank,
         max_hu_bank)


sa_summary = sa_summary %>% 
  dplyr::rename(id = sa_id,
                sum_hu_sa = sum_housing_units,
                mean_hu_sa = mean_housing_units,
                median_hu_sa = median_housing_units,
                max_hu_sa = max_housing_units) %>% 
  dplyr::select(id, 
                sum_hu_sa,
                mean_hu_sa,
                median_hu_sa,
                max_hu_sa)

loss_summary = loss_summary %>% 
  dplyr::rename(id = sa_id,
                sum_hu_loss = sum_housing_units,
                mean_hu_loss = mean_housing_units,
                median_hu_loss = median_housing_units,
                max_hu_loss = max_housing_units) %>% 
  dplyr::select(id, 
                sum_hu_loss,
                mean_hu_loss,
                median_hu_loss,
                max_hu_loss)

summary <- bank_summary %>%
  inner_join(sa_summary, by = "id") %>%
  inner_join(loss_summary, by = "id")
