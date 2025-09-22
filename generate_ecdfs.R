# Packages ----------------------------------------------------
library(future)
library(future.apply)
library(dplyr)
library(purrr)
library(tidyr)

# Base directory (override with env‑var) ----------------------
datadir <- Sys.getenv("DATA_DIR", "C:/Users/indumati/Box/Paper2_final")

# Helper to build absolute paths quickly
pth <- function(...) file.path(datadir, ...)

#-----------------------------------------------------#
#     Saving ECDF Functions for 2021 Service Areas    #
#-----------------------------------------------------#

extracted_dir <- pth("Extractions and Summaries", "Extract Service Areas")
output_dir_sa <- pth("Extractions and Summaries", "ECDF Functions", "SA ECDF Functions")
dir.create(output_dir_sa, recursive = TRUE, showWarnings = FALSE)

rds_files <- list.files(extracted_dir, pattern = "^sa_extract_.*\\.rds$", full.names = TRUE)

existing_ecdfs <- list.files(output_dir_sa, pattern = "^ecdf_fn_.*\\.rds$", full.names = FALSE)
processed_sa_names <- gsub("^ecdf_fn_(.*)\\.rds$", "\\1", existing_ecdfs)
unprocessed_files <- rds_files[!gsub("^sa_extract_(.*)\\.rds$", "\\1", basename(rds_files)) %in% processed_sa_names]

compute_ecdf_functions <- function(df, sa_name) {
  bank_extract <- readRDS(pth("Extractions and Summaries", "Extract Banks", "bank_extractions.rds"))
  bank_extract_sa <- bank_extract[bank_extract$bank_name == sa_name, ]
  df <- df[!df$cell %in% bank_extract_sa$cell, ]
  
  vars <- c("housing_units", "housing_value", "population")
  ecdf_func_list <- lapply(vars, function(var) {
    if (!var %in% names(df)) return(NULL)
    valid_values <- na.omit(df[[var]])
    if (length(valid_values) == 0) return(NULL)
    ecdf(valid_values)
  })
  names(ecdf_func_list) <- vars
  return(ecdf_func_list)
}

plan(multisession, workers=10)

process_sa <- function(file_path) {
  sa_name <- gsub("^sa_extract_(.*)\\.rds$", "\\1", basename(file_path))
  if (sa_name %in% processed_sa_names) return(NULL)
  df <- tryCatch(readRDS(file_path), error = function(e) return(NULL))
  if (is.null(df)) return(NULL)
  if (ncol(df) >= 4) colnames(df)[2:4] <- c("housing_units", "housing_value", "population")
  ecdf_funcs <- compute_ecdf_functions(df, sa_name)
  saveRDS(ecdf_funcs, file.path(output_dir_sa, paste0("ecdf_fn_", sa_name, ".rds")))
  return(sa_name)
}

results <- future_lapply(unprocessed_files, process_sa)

#-----------------------------------------------------#
#        Saving ECDF Functions for 2021 Banks         #
#-----------------------------------------------------#

extracted_dir <- pth("Extractions and Summaries", "Extract Banks")
output_dir_banks <- pth("Extractions and Summaries", "ECDF Functions", "Bank ECDF Functions")
dir.create(output_dir_banks, recursive = TRUE, showWarnings = FALSE)

rds_files <- list.files(extracted_dir, pattern = "^bank_extract_.*\\.rds$", full.names = TRUE)
existing_ecdfs <- list.files(output_dir_banks, pattern = "^bank_ecdf_fn_.*\\.rds$", full.names = FALSE)
processed_bank_names <- gsub("^bank_ecdf_fn_(.*)\\.rds$", "\\1", existing_ecdfs)
unprocessed_files <- rds_files[!gsub("^bank_extract_(.*)\\.rds$", "\\1", basename(rds_files)) %in% processed_bank_names]

compute_ecdf_functions <- function(df) {
  vars <- c("housing_units", "housing_value", "population")
  ecdf_func_list <- lapply(vars, function(var) {
    if (!var %in% names(df)) return(NULL)
    valid_values <- na.omit(df[[var]])
    if (length(valid_values) == 0) return(NULL)
    ecdf(valid_values)
  })
  names(ecdf_func_list) <- vars
  return(ecdf_func_list)
}

process_bank <- function(file_path) {
  bank_name <- gsub("^bank_extract_(.*)\\.rds$", "\\1", basename(file_path))
  if (bank_name %in% processed_bank_names) return(NULL)
  df <- tryCatch(readRDS(file_path), error = function(e) return(NULL))
  if (is.null(df)) return(NULL)
  if (ncol(df) >= 4) colnames(df)[2:4] <- c("housing_units", "housing_value", "population")
  ecdf_funcs <- compute_ecdf_functions(df)
  saveRDS(ecdf_funcs, file.path(output_dir_banks, paste0("bank_ecdf_fn_", bank_name, ".rds")))
  return(bank_name)
}

plan(multisession, workers = 4)
results <- future_lapply(unprocessed_files, process_bank)

#-----------------------------------------------------#
#   Saving ECDF Functions for Service Areas - Loss    #
#-----------------------------------------------------#

extracted_dir <- pth("Extractions and Summaries", "Extract Loss")
output_dir_sa <- pth("Extractions and Summaries", "ECDF Functions", "Loss ECDF Functions")
dir.create(output_dir_sa, recursive = TRUE, showWarnings = FALSE)

rds_files <- list.files(extracted_dir, pattern = "^loss_extract_.*\\.rds$", full.names = TRUE)
existing_ecdfs <- list.files(output_dir_sa, pattern = "^ecdf_fn_.*\\.rds$", full.names = FALSE)
processed_sa_names <- gsub("^ecdf_fn_(.*)\\.rds$", "\\1", existing_ecdfs)
unprocessed_files <- rds_files[!gsub("^loss_extract_(.*)\\.rds$", "\\1", basename(rds_files)) %in% processed_sa_names]

compute_ecdf_functions <- function(df, sa_name) {
  bank_extract <- readRDS(pth("Extractions and Summaries", "Extract Bank", "bank_extractions.rds"))
  df <- df[!df$cell %in% bank_extract[bank_extract$bank_name == sa_name, "cell"], ]
  vars <- c("housing_units", "housing_value", "population")
  ecdf_func_list <- lapply(vars, function(var) {
    if (!var %in% names(df)) return(NULL)
    valid_values <- na.omit(df[[var]])
    if (length(valid_values) == 0) return(NULL)
    ecdf(valid_values)
  })
  names(ecdf_func_list) <- vars
  return(ecdf_func_list)
}

plan(multisession, workers = 10)
process_sa <- function(file_path) {
  sa_name <- gsub("^sa_extract_(.*)\\.rds$", "\\1", basename(file_path))
  if (sa_name %in% processed_sa_names) return(NULL)
  df <- tryCatch(readRDS(file_path), error = function(e) return(NULL))
  if (is.null(df)) return(NULL)
  if (ncol(df) >= 4) colnames(df)[2:4] <- c("housing_units", "housing_value", "population")
  ecdf_funcs <- compute_ecdf_functions(df, sa_name)
  saveRDS(ecdf_funcs, file.path(output_dir_sa, paste0("ecdf_fn_", sa_name, ".rds")))
  return(sa_name)
}

results <- future_lapply(unprocessed_files, process_sa)

#-----------------------------------------------------#
#          Combine All ECDFs into Final Table         #
#-----------------------------------------------------#

bank2021_ecdf_dir <- pth("Extractions and Summaries", "ECDF Functions", "Bank ECDF Functions")
sa2021_ecdf_dir <- pth("Extractions and Summaries", "ECDF Functions", "SA ECDF Functions")
sa_loss_ecdf_dir <- pth("Extractions and Summaries", "ECDF Functions", "Loss ECDF Functions")

bank2021_files <- list.files(bank2021_ecdf_dir, pattern = "^bank_ecdf_fn_.*\\.rds$", full.names = TRUE)
sa2021_files <- list.files(sa2021_ecdf_dir, pattern = "^ecdf_fn_.*\\.rds$", full.names = TRUE)
sa_loss_files <- list.files(sa_loss_ecdf_dir, pattern = "^ecdf_fn_.*\\.rds$", full.names = TRUE)

parse_bank_name <- function(filename) sub("^bank_ecdf_fn_(.*)\\.rds$", "\\1", basename(filename))
parse_sa_name <- function(filename) sub("^ecdf_fn_(.*)\\.rds$", "\\1", basename(filename))
read_ecdf_list <- function(file_path) readRDS(file_path)

make_ecdf_tibble <- function(files, name_parser, tag) {
  tibble(
    file = files,
    entity = map_chr(files, name_parser),
    ecdf_list = map(files, read_ecdf_list)
  ) %>%
    mutate(
      housing_units = map(ecdf_list, ~ .x[["housing_units"]] %||% NA),
      housing_value = map(ecdf_list, ~ .x[["housing_value"]] %||% NA),
      population = map(ecdf_list, ~ .x[["population"]] %||% NA)
    ) %>%
    select(entity, housing_units, housing_value, population) %>%
    rename_with(~ paste0(tag, "_", .), -entity)
}

bank2021_tbl <- make_ecdf_tibble(bank2021_files, parse_bank_name, "Bank_2021")
sa2021_tbl <- make_ecdf_tibble(sa2021_files, parse_sa_name, "SA_2021")
sa_loss_tbl <- make_ecdf_tibble(sa_loss_files, parse_sa_name, "SA_loss")

common_entities <- Reduce(intersect, list(
  bank2021_tbl$entity,
  sa2021_tbl$entity,
  sa_loss_tbl$entity
))

bank2021_tbl <- filter(bank2021_tbl, entity %in% common_entities)
sa2021_tbl <- filter(sa2021_tbl, entity %in% common_entities)
sa_loss_tbl <- filter(sa_loss_tbl, entity %in% common_entities)

final_tbl <- bank2021_tbl %>%
  full_join(sa2021_tbl, by = "entity") %>%
  full_join(sa_loss_tbl, by = "entity")

final_tbl_long <- final_tbl %>%
  pivot_longer(-entity, names_to = c("data_type", "variable"),
               names_pattern = "^(.*)_(housing_units|housing_value|population)$",
               values_to = "ecdf_func") %>%
  pivot_wider(names_from = "variable", values_from = "ecdf_func", names_glue = "{variable}_ecdf") %>%
  arrange(entity, data_type)
