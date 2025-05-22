library(future)
library(future.apply)
library(dplyr)


#-----------------------------------------------------#
#     Saving ECDF Functions for 2021 Service Areas    #
#-----------------------------------------------------#

# Setwd(your_directory)

# Get names of service areas from the extraction files
extracted_dir = "Summaries and Extractions/ Extract Service Areas"
rds_files <- list.files(
extracted_dir, 
pattern = "^sa_extract_.*\\.rds$", 
full.names = TRUE
)


output_dir_sa = "ECDF Functions/SA ECDF Functions"


existing_ecdfs <- list.files(
  output_dir_sa, 
  pattern = "^ecdf_fn_.*\\.rds$", 
  full.names = FALSE
)

processed_sa_names <- gsub("^ecdf_fn_(.*)\\.rds$", "\\1", existing_ecdfs)
unprocessed_files <- rds_files[!gsub("^sa_extract_(.*)\\.rds$", "\\1", basename(rds_files)) %in% processed_sa_names]


compute_ecdf_functions <- function(df, sa_name) {
  
  # Read bank extract
  bank_extract <- readRDS("Summaries and Extractions/Extract Banks/all_bank_extractions.rds")
  
  # Subset to only the banks associated with this service area
  bank_extract_sa <- bank_extract[bank_extract$bank_name == sa_name, ]
  
  # Extract bank cell IDs
  bank_cell_ids <- bank_extract_sa$cell
  
  # Remove cells that appear in bank_extract
  df <- df[!df$cell %in% bank_cell_ids, ]
  
  vars <- c("housing_units", "housing_value", "population")
  
  ecdf_func_list <- list()
  
  for (var in vars) {
    if (!var %in% colnames(df)) next
    valid_values <- na.omit(df[[var]])
    if (length(valid_values) == 0) {
      ecdf_func_list[[var]] <- NULL
    } else {
      ecdf_func_list[[var]] <- ecdf(valid_values)
    }
  }
  
  return(ecdf_func_list)
}


plan(multisession, workers=10)  # Parallelize

process_sa <- function(file_path) {
  
  base_name <- basename(file_path)
  sa_name   <- gsub("^sa_extract_(.*)\\.rds$", "\\1", base_name)
  
  if (sa_name %in% processed_sa_names) {
    cat("\nSkipping already processed:", sa_name, "\n")
    return(NULL)
  }
  
  cat("\n---------------------------------------------------\n")
  cat("Processing Service Area:", sa_name, "\n")
  
  df <- tryCatch({
    readRDS(file_path)
  }, error = function(e) {
    cat("  ERROR reading file:", file_path, "\n  Message:", e$message, "\n  Skipping.\n")
    return(NULL)
  })
  
  if (is.null(df)) return(NULL)
  
  if (ncol(df) >= 4) {
    colnames(df)[2:4] <- c("housing_units", "housing_value", "population")
  }
  
  # Compute ECDF functions
  ecdf_funcs  <- compute_ecdf_functions(df, sa_name)
  
  # Save ECDF function objects
  ecdf_rds_path <- file.path(output_dir_sa, paste0("ecdf_fn_", sa_name, ".rds"))
  saveRDS(ecdf_funcs, ecdf_rds_path)
  
  cat("  ECDF functions saved:", ecdf_rds_path, "\n")
  
  return(sa_name)
}

# Process files in parallel
results <- future_lapply(unprocessed_files, process_sa)



#-----------------------------------------------------#
#        Saving ECDF Functions for 2021 Banks         #
#-----------------------------------------------------#


# Where the bank_extract_*.rds files live
extracted_dir <- "Extractions and Summaries/Extract Banks"

# Where to store the final ECDF function outputs
output_dir_banks <- file.path("ECDF Functions/Bank ECDF Functions")
dir.create(output_dir_banks, recursive = TRUE, showWarnings = FALSE)


# Find all bank_extract_*.rds, and which are processed

rds_files <- list.files(
  extracted_dir, 
  pattern = "^bank_extract_.*\\.rds$", 
  full.names = TRUE
)

existing_ecdfs <- list.files(
  output_dir_banks, 
  pattern = "^bank_ecdf_fn_.*\\.rds$", 
  full.names = FALSE
)

processed_bank_names <- gsub("^bank_ecdf_fn_(.*)\\.rds$", "\\1", existing_ecdfs)

# filter out any bank_extract files for already processed banks
unprocessed_files <- rds_files[ 
  !gsub("^bank_extract_(.*)\\.rds$", "\\1", basename(rds_files)) %in% processed_bank_names
]


#  function to compute ecdf() across 3 columns

compute_ecdf_functions <- function(df) {
  vars <- c("housing_units", "housing_value", "population")
  ecdf_func_list <- list()
  
  for (var in vars) {
    if (!var %in% colnames(df)) next
    
    valid_values <- na.omit(df[[var]])
    if (length(valid_values) == 0) {
      ecdf_func_list[[var]] <- NULL
    } else {
      ecdf_func_list[[var]] <- ecdf(valid_values)
    }
  }
  
  return(ecdf_func_list)
}

# Processing function for a single bank_extract file

process_bank <- function(file_path) {
  
  base_name <- basename(file_path)
  bank_name <- gsub("^bank_extract_(.*)\\.rds$", "\\1", base_name)
  
  # skip if bank_name already processed
  if (bank_name %in% processed_bank_names) {
    cat("\nSkipping already processed:", bank_name, "\n")
    return(NULL)
  }
  
  cat("\n---------------------------------------------------\n")
  cat("Processing Bank:", bank_name, "\n")
  
  # Read the extracted data
  df <- tryCatch({
    readRDS(file_path)
  }, error = function(e) {
    cat("  ERROR reading file:", file_path, "\n  Message:", e$message, "\n  Skipping.\n")
    return(NULL)
  })
  
  # If reading failed, skip
  if (is.null(df)) return(NULL)
  
  # Rename columns 2,3,4 -> housing_units, housing_value, population
  
  if (ncol(df) >= 4) {
    colnames(df)[2:4] <- c("housing_units", "housing_value", "population")
  }
  
  # Compute ECDF function objects
  ecdf_funcs <- compute_ecdf_functions(df)
  
  # Save the result
  ecdf_rds_path <- file.path(output_dir_banks, paste0("bank_ecdf_fn_", bank_name, ".rds"))
  saveRDS(ecdf_funcs, ecdf_rds_path)
  
  cat("  ECDF functions saved:", ecdf_rds_path, "\n")
  
  return(bank_name)
}


# Parallelize and run

plan(multisession, workers = 2)

results <- future_lapply(unprocessed_files, process_bank)

cat("\nDone! Processed", length(results), "banks.\n")



#-----------------------------------------------------#
#        Saving ECDF Functions for 2021 Banks         #
#-----------------------------------------------------#

extracted_dir <- "CONUS/Wetland Loss/Extractions/Extract Service Area Loss"             # Where the .rds extracted data live
output_dir_sa <- file.path(extracted_dir, "SA Summaries Loss/SA ECDF Functions")            # Where to save the summary files
dir.create(output_dir_sa, recursive = TRUE, showWarnings = FALSE)

rds_files <- list.files(
  extracted_dir, 
  pattern = "^sa_extract_.*\\.rds$", 
  full.names = TRUE
)
existing_ecdfs <- list.files(
  output_dir_sa, 
  pattern = "^ecdf_fn_.*\\.rds$", 
  full.names = FALSE
)

processed_sa_names <- gsub("^ecdf_fn_(.*)\\.rds$", "\\1", existing_ecdfs)
unprocessed_files <- rds_files[!gsub("^sa_extract_(.*)\\.rds$", "\\1", basename(rds_files)) %in% processed_sa_names]


compute_ecdf_functions <- function(df, sa_name) {
  
  # Read bank extract
  bank_extract <- readRDS("CONUS/Wetland Loss/Extractions/Extract Bank Loss/bank_extractions_wetlandloss.rds")
  
  # Subset to only the banks associated with this service area
  bank_extract_sa <- bank_extract[bank_extract$bank_name == sa_name, ]
  
  # Extract bank cell IDs
  bank_cell_ids <- bank_extract_sa$cell
  
  # Remove cells in bank_extract
  df <- df[!df$cell %in% bank_cell_ids, ]
  
  vars <- c("housing_units", "housing_value", "population")
  
  ecdf_func_list <- list()
  
  for (var in vars) {
    if (!var %in% colnames(df)) next
    valid_values <- na.omit(df[[var]])
    if (length(valid_values) == 0) {
      ecdf_func_list[[var]] <- NULL
    } else {
      ecdf_func_list[[var]] <- ecdf(valid_values)
    }
  }
  
  return(ecdf_func_list)
}

# parallelize

plan(multisession, workers=10)

process_sa <- function(file_path) {
  
  base_name <- basename(file_path)
  sa_name   <- gsub("^sa_extract_(.*)\\.rds$", "\\1", base_name)
  
  if (sa_name %in% processed_sa_names) {
    cat("\nSkipping already processed:", sa_name, "\n")
    return(NULL)
  }
  
  cat("\n---------------------------------------------------\n")
  cat("Processing Service Area:", sa_name, "\n")
  
  df <- tryCatch({
    readRDS(file_path)
  }, error = function(e) {
    cat("  ERROR reading file:", file_path, "\n  Message:", e$message, "\n  Skipping.\n")
    return(NULL)
  })
  
  if (is.null(df)) return(NULL)
  
  if (ncol(df) >= 4) {
    colnames(df)[2:4] <- c("housing_units", "housing_value", "population")
  }
  
  # Compute ECDF functions
  ecdf_funcs  <- compute_ecdf_functions(df, sa_name)
  
  # Save ECDF function objects
  ecdf_rds_path <- file.path(output_dir_sa, paste0("ecdf_fn_", sa_name, ".rds"))
  saveRDS(ecdf_funcs, ecdf_rds_path)
  
  cat("  ECDF functions saved:", ecdf_rds_path, "\n")
  
  return(sa_name)
}

results <- future_lapply(unprocessed_files, process_sa)

###################################################################
#                                                                 #
#               Combining into one dataset                        #
#                                                                 #
###################################################################

library(dplyr)
library(purrr)
library(tidyr)


#-----------------------------------------------------------------------------
#  Get directories where the ECDF .rds files live (adjust ofc)
#-----------------------------------------------------------------------------
bank2021_ecdf_dir  <- "Extractions/Extract Banks 2021/Bank ECDF Functions"

sa2021_ecdf_dir    <- "Extractions/Extract Service Areas 2021/SA ECDF Functions"  

sa_loss_ecdf_dir   <- "CONUS/Wetland Loss/Extractions/Extract Service Area Loss/SA Summaries Loss/SA ECDF Functions"

#-----------------------------------------------------------------------------
#  Get file paths
#-----------------------------------------------------------------------------
bank2021_files <- list.files(
  bank2021_ecdf_dir, 
  pattern     = "^bank_ecdf_fn_.*\\.rds$", 
  full.names  = TRUE
)

sa2021_files <- list.files(
  sa2021_ecdf_dir,
  pattern     = "^ecdf_fn_.*\\.rds$",
  full.names  = TRUE
)

sa_loss_files <- list.files(
  sa_loss_ecdf_dir,
  pattern     = "^ecdf_fn_.*\\.rds$",
  full.names  = TRUE
)

#-----------------------------------------------------------------------------
# Define Helper functions
#-----------------------------------------------------------------------------
parse_bank_name <- function(filename) {
  # bank_ecdf_fn_<BANKNAME>.rds
  sub("^bank_ecdf_fn_(.*)\\.rds$", "\\1", basename(filename))
}

parse_sa_name <- function(filename) {
  # ecdf_fn_<SANAME>.rds
  sub("^ecdf_fn_(.*)\\.rds$", "\\1", basename(filename))
}

read_ecdf_list <- function(file_path) {
  # Each RDS is a list of 3 ecdf objects:
  #   $housing_units, $housing_value, $population
  readRDS(file_path)
}

#-----------------------------------------------------------------------------
#  Create tibble for each dataset type
#-----------------------------------------------------------------------------
make_ecdf_tibble <- function(files, name_parser, tag) {
    tbl <- tibble(
    file   = files,
    entity = map_chr(files, name_parser),
    ecdf_list = map(files, read_ecdf_list)
  )
  
  #  Extract each of the 3 ecdf objects into a separate list-column
  tbl <- tbl %>%
    mutate(
      housing_units = map(ecdf_list, ~ .x[["housing_units"]] %||% NA),
      housing_value = map(ecdf_list, ~ .x[["housing_value"]] %||% NA),
      population    = map(ecdf_list, ~ .x[["population"]]    %||% NA)
    ) %>%
    select(entity, housing_units, housing_value, population)
  
  # 3) Rename as wide format
  tbl <- tbl %>%
    rename_with(~ paste0(tag, "_", .), -entity)
  
  return(tbl)
}

bank2021_tbl  <- make_ecdf_tibble(bank2021_files, parse_bank_name, "Bank_2021")
sa2021_tbl    <- make_ecdf_tibble(sa2021_files,   parse_sa_name,  "SA_2021")
sa_loss_tbl   <- make_ecdf_tibble(sa_loss_files,  parse_sa_name,  "SA_loss")

# 
bad_files <- c()
for (f in sa2021_files) {
  rds_data <- tryCatch({
    readRDS(f)
  },
  error = function(e) {
    bad_files <<- c(bad_files, f)
    NULL
  }
  )
}
bad_files

final_tbl <- bank2021_tbl %>%
  full_join(sa2021_tbl, by = "entity") %>%
  full_join(sa_loss_tbl, by = "entity")



########### GET ONLY THOSE IN ALL DATASETS:
# Get entity names from each
bank_entities    <- bank2021_tbl$entity
sa2021_entities  <- sa2021_tbl$entity
sa_loss_entities <- sa_loss_tbl$entity

common_entities <- Reduce(intersect, list(bank_entities, sa2021_entities, sa_loss_entities))

length(common_entities)  

bank2021_tbl <- bank2021_tbl %>%
  filter(entity %in% common_entities)

sa2021_tbl <- sa2021_tbl %>%
  filter(entity %in% common_entities)

sa_loss_tbl <- sa_loss_tbl %>%
  filter(entity %in% common_entities)

final_tbl <- bank2021_tbl %>%
  full_join(sa2021_tbl, by = "entity") %>%
  full_join(sa_loss_tbl, by = "entity")


# Pivot longer 
final_tbl_long <- final_tbl %>%
  pivot_longer(
    cols = -entity,
    names_to = c("data_type", "variable"),
    names_pattern = "^(.*)_(housing_units|housing_value|population)$",
    values_to = "ecdf_func"
  ) %>%
  pivot_wider(
    names_from = "variable",
    values_from = "ecdf_func",
    names_glue = "{variable}_ecdf"
  ) %>%
  arrange(entity, data_type)



