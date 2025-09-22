suppressPackageStartupMessages({
  library(tidyverse)
  library(ggplot2)
  library(scales)
  library(patchwork)
  library(vctrs)
  library(xml2)
  library(ggplot2)
  library(terra)
  library(sf)
  library(tmap)
  library(tidycensus)
  library(data.table)
  library(ggspatial)
  library(future)
  library(future.apply)
  library(dplyr)
  library(purrr)
  library(tidyr)
  library(fst)
  library(progressr)
  library(reasdr)
})


pth <- function(...) file.path("C:/Users/indumati/Box/Paper2_final", ...)

#-----------------------------------------------------#
# Load banks and service areas
#-----------------------------------------------------#
setwd("C:/Users/indumati/Box/Paper2_final")

serviceareas = readRDS(pth("Service Areas", "ServiceAreas_agg.rds"))
banks=readRDS(pth("Footprints", "footprints_and_buffers.rds"))

banks <- st_as_sf(banks)
sas <- st_as_sf(sas)

banks = banks %>%
  mutate(Name = recode(Name, 
                       "FP_AndL_Everglades_Phase_II_MB" = "FP_L_Everglades_Phase_II_MB"))

sanitize_name <- function(name) {
  name <- gsub("[^A-Za-z0-9_-]", "_", name)
  name <- gsub("_+", "_", name)
  name <- gsub("^_|_$", "", name)
  return(name)
}

summary_bank = readRDS(pth("Extractions and Summaries", "bank_summary.rds"))

banks = banks %>% 
  select(c("Name", "geometry"))
sas <- sas %>% 
  select(c("ID", "geometry"))

sas <- sas %>%
  filter(ID %in% summary_bank$bank_name)
banks <- banks %>%
  filter(Name %in% summary_bank$bank_name)


# FIND IDENTICAL SAS
library(purrr)
sas_bbox <- sas %>%
  st_make_valid() %>%
  mutate(
    bbox = map(geometry, ~st_bbox(.)),
    xmin = map_dbl(bbox, ~.[1]),
    ymin = map_dbl(bbox, ~.[2]),
    xmax = map_dbl(bbox, ~.[3]),
    ymax = map_dbl(bbox, ~.[4]),
    # Round to nearest 10 units instead of 0.1
    xmin_r = round(xmin, -1),  # -1 rounds to nearest 10
    ymin_r = round(ymin, -1),
    xmax_r = round(xmax, -1),
    ymax_r = round(ymax, -1)
  )

bbox_duplicates <- sas_bbox %>%
  st_drop_geometry() %>%
  group_by(xmin_r, ymin_r, xmax_r, ymax_r) %>%
  filter(n() > 1) %>%
  summarise(
    identical_sa_ids = paste(ID, collapse = ", "),
    count = n(),
    .groups = "drop"
  )

print(bbox_duplicates)


#-----------------------------------------------------#
# Pooled Bank ECDF Functions for Identical Service Areas
#-----------------------------------------------------#

########

# ------------------------------------------------------------
#               Libraries & helper utilities
# ------------------------------------------------------------
sanitize_name <- function(x) {
  x %>% 
    str_replace_all("[^A-Za-z0-9_-]", "_") %>%   # keep safe chars
    str_replace_all("_+", "_")  %>%              # collapse repeats
    str_remove_all("^_|_$")                      # trim ends
}

compute_ecdf_set <- function(df, vars = c("housing_units",
                                          "housing_value",
                                          "population")) {
  out <- list()
  for (v in vars) {
    if (!v %in% names(df)) next
    vals <- na.omit(df[[v]])
    out[[v]] <- if (length(vals)) ecdf(vals) else NULL
  }
  out
}

# ------------------------------------------------------------
#            1.  pool ⇄ bank cross-walk  (once)
# ------------------------------------------------------------
bank_pool_map <- bbox_duplicates %>% 
  mutate(pool_id = paste0("Pool_", row_number())) %>% 
  select(pool_id, identical_sa_ids) %>% 
  separate_rows(identical_sa_ids, sep = ",\\s*") %>% 
  rename(bank_name_raw = identical_sa_ids) %>% 
  mutate(bank_name = sanitize_name(bank_name_raw))

# ------------------------------------------------------------
#            2.  locate extract files on disk 
# ------------------------------------------------------------

bank_extract_dir <- pth("CONUS", "Extractions and Summaries", "Extract Banks")
sa_extract_dir   <- pth("CONUS","Extractions and Summaries", "Extract Service Areas")

bank_files <- tibble(file_path = list.files(bank_extract_dir,
                                            pattern = "^bank_extract_.*\\.rds$",
                                            full.names = TRUE)) %>% 
  mutate(bank_name = gsub("^bank_extract_(.*)\\.rds$", "\\1",
                          basename(file_path))) %>% 
  left_join(bank_pool_map, by = "bank_name")

sa_files <- tibble(file_path = list.files(sa_extract_dir,
                                          pattern = "^sa_extract_.*\\.rds$",
                                          full.names = TRUE)) %>% 
  mutate(sa_id = gsub("^sa_extract_(.*)\\.rds$", "\\1",
                      basename(file_path))) %>% 
  left_join(bank_pool_map %>% 
              distinct(pool_id, sa_id = bank_name),   # same cleaned label
            by = "sa_id")

# ------------------------------------------------------------
#            3.  main worker: one pool at a time
# ------------------------------------------------------------
rename_layers <- function(df) {
  if (ncol(df) >= 4)               
    colnames(df)[2:4] <- c("housing_units",
                           "housing_value",
                           "population")
  df
}

compare_one_pool <- function(pid,
                             vars = c("housing_units",
                                      "housing_value",
                                      "population")) {
  
  sa_path    <- sa_files  %>% filter(pool_id == pid) %>% pull(file_path) %>% first()
  bank_paths <- bank_files %>% filter(pool_id == pid) %>% pull(file_path)
  
  if (length(sa_path) == 0 || length(bank_paths) == 0)
    stop("No SA or bank extracts found for pool ", pid, call. = FALSE)
  
  # --- load & rename --------------------------------------------------------
  sa_dat   <- rename_layers(read_rds(sa_path))
  bank_dfs <- map(bank_paths, ~ rename_layers(read_rds(.x)))
  bank_dat <- bind_rows(bank_dfs)
  
  # ---- mask out bank cells from SA ----------------------------------------
  bank_cells <- unique(bank_dat$cell)
  sa_remainder <- sa_dat[ !(sa_dat$cell %in% bank_cells), ]
  
  # ---- ECDFs ---------------------------------------------------------------
  ecdf_bank   <- compute_ecdf_set(bank_dat, vars)
  ecdf_sa_rem <- compute_ecdf_set(sa_remainder, vars)
  
  list(
    pool_id        = pid,
    n_bank_cells   = nrow(bank_dat),
    n_sa_cells     = nrow(sa_dat),
    n_sa_rem_cells = nrow(sa_remainder),
    bank_ecdf      = ecdf_bank,
    sa_rem_ecdf    = ecdf_sa_rem,
    meta = list(
      bank_files = bank_paths,
      sa_file    = sa_path,
      processed  = Sys.time()
    )
  )
}


# ------------------------------------------------------------
#            4.  run for every pool & save
# ------------------------------------------------------------
all_pools   <- sort(unique(bank_pool_map$pool_id))
output_dir  <- pth("ECDF Functions", "Pool_vs_SA_remainder")
dir.create(output_dir, recursive = TRUE, showWarnings = FALSE)

results <- map(all_pools, \(pid) {
  message("Processing ", pid, " …")
  res <- compare_one_pool(pid)
  write_rds(res, file.path(output_dir,
                           paste0("ecdf_pool_vs_sa_", pid, ".rds")))
  invisible(res)   # keeps memory light
})

message("Done!  Results saved in: ", output_dir)

#----------------------------------
# %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

# Make PDF plot


library(tidyverse)
library(scales)
library(vctrs)
library(tidyr)
library(stringr)

compare_remaining_vs_lost <- function(datadir, sample_size = 10000, seed = 123, bank_filter = NULL) {
  
  message("Scanning directory for bank ECDF files …")
  bank_files <- list.files(file.path(datadir, "SA ECDF Functions"), full.names = FALSE)
  bank_names <- sub("^ecdf_fn_(.*)\\.rds$", "\\1", bank_files)
  
  # ✅ Filter bank_names if filter is provided
  if (!is.null(bank_filter)) {
    bank_names <- bank_names[bank_names %in% bank_filter]
  }
  
  total_banks <- length(bank_names)

  ## --- Service area duplicate weighting ---
  sa_lookup <- bbox_duplicates %>%
    separate_rows(identical_sa_ids, sep = ",\\s*") %>%
    mutate(bank = str_trim(identical_sa_ids)) %>%
    select(bank, count)
  
  count_lookup <- tibble(bank = bank_names) %>%
    left_join(sa_lookup, by = "bank") %>%
    mutate(count = replace_na(count, 1L))
  
  ## --- Set up ---
  results_list <- vector("list", total_banks * 2)
  j <- 1
  successful_banks <- 0
  set.seed(seed)
  
  for (i in seq_along(bank_names)) {
    bank <- bank_names[i]
    message(sprintf("[%d/%d] Processing bank: %s", i, total_banks, bank))
    
    sa_file   <- file.path(datadir, "SA ECDF Functions",        sprintf("ecdf_fn_%s.rds", bank))
    loss_file <- file.path(datadir, "SA ECDF Functions - Loss", sprintf("ecdf_fn_%s.rds", bank))
    
    if (!file.exists(sa_file) || !file.exists(loss_file)) {
      message("  - Skipping – missing ECDF file(s)")
      next
    }
    
    sa_ecdf   <- readRDS(sa_file)$housing_units
    loss_ecdf <- readRDS(loss_file)$housing_units
    if (is.null(sa_ecdf) || is.null(loss_ecdf)) {
      message("  - Skipping – housing_units ECDF missing")
      next
    }
    
    sa_data   <- knots(sa_ecdf)
    loss_data <- knots(loss_ecdf)
    if (length(sa_data) == 0 || length(loss_data) == 0) {
      message("  - Skipping – ECDF contains no data")
      next
    }
    
    sa_n_obs   <- environment(sa_ecdf)$n   %||% length(sa_data)
    loss_n_obs <- environment(loss_ecdf)$n %||% length(loss_data)
    
    # Lookup how many banks share this service area
    dup_factor <- count_lookup$count[count_lookup$bank == bank]
    row_weight <- (loss_n_obs / length(loss_data)) / dup_factor
    
    bank_remaining_mean <- mean(sa_data)
    
    if (length(sa_data) > sample_size) {
      sa_data <- sample(sa_data, sample_size)
      message(sprintf("  - Downsampled remaining: %s -> %s",
                      comma(sa_n_obs), comma(length(sa_data))))
    }
    if (length(loss_data) > sample_size) {
      loss_data <- sample(loss_data, sample_size)
      message(sprintf("  - Downsampled lost: %s -> %s",
                      comma(loss_n_obs), comma(length(loss_data))))
    }
    
    remaining <- tibble(
      bank, status = "remaining", value = sa_data,
      normalized_value = sa_data / bank_remaining_mean,
      weight = row_weight
    ) %>% filter(is.finite(normalized_value) & normalized_value > 0)
    
    lost <- tibble(
      bank, status = "lost", value = loss_data,
      normalized_value = loss_data / bank_remaining_mean,
      weight = row_weight
    ) %>% filter(is.finite(normalized_value) & normalized_value > 0)
    
    results_list[[j]] <- remaining; j <- j + 1
    results_list[[j]] <- lost;      j <- j + 1
    successful_banks  <- successful_banks + 1
  }
  
  if (successful_banks == 0) stop("No banks with both remaining and lost ECDFs were processed.")
  message("Successfully processed ", successful_banks, " bank(s)")
  
  results <- bind_rows(results_list)
  rm(results_list); gc()
  
  message("Computing ECDF data …")
  results_ecdf_full <- results %>%
    group_by(bank, status) %>%
    arrange(normalized_value) %>%
    mutate(ecdf_y = row_number() / n()) %>%
    ungroup()
  
  ## ---------- statistics for the PDF panel ----------
  loss_mean <- weighted.mean(results$normalized_value[results$status == "lost"],
                             w = results$weight[results$status == "lost"], na.rm = TRUE)
  
  # t-test (unweighted); 
  p_val <- t.test(normalized_value ~ status, data = results)$p.value
  sig_stars <- case_when(
    p_val < 0.001 ~ "***",
    p_val < 0.01  ~ "**",
    p_val < 0.05  ~ "*",
    TRUE          ~ "ns"
  )
  
  # Rough y-coordinate for bracket
  dens_rem  <- density(results$normalized_value[results$status == "remaining"])
  dens_lost <- density(results$normalized_value[results$status == "lost"])
  y_max     <- max(dens_rem$y, dens_lost$y, na.rm = TRUE)
  y_bracket <- y_max * 1.05
  y_tick    <- y_max * 0.97
  x_left  <- min(1, loss_mean)
  x_right <- max(1, loss_mean)
  
  log_breaks <- c(0.1, 1, 10, 100)
  
  pdf_plot <- ggplot(results_ecdf_full,
                     aes(x = normalized_value,
                         colour = status,
                         weight = weight)) +
    stat_density(geom = "line", position = "identity", linewidth = 1.2) +
    geom_vline(xintercept = 1,         linetype = "dashed", colour = "blue", linewidth = .8) +
    geom_vline(xintercept = loss_mean, linetype = "dashed", colour = "red",  linewidth = .8) +
    annotate("segment", x = x_left,  xend = x_right, y = y_bracket, yend = y_bracket) +
    annotate("segment", x = x_left,  xend = x_left,  y = y_bracket, yend = y_tick) +
    annotate("segment", x = x_right, xend = x_right, y = y_bracket, yend = y_tick) +
    annotate("text", x = (x_left + x_right)/2 - 0.2, y = y_bracket + 0.03 * y_max,
             label = sig_stars, vjust = 0)+
    scale_x_log10(
      limits = c(0.01, 100),
      breaks = sort(unique(c(0.1, 1, 10, 100, signif(loss_mean, 2)))),
      labels = function(x) format(x, digits = 2, trim = TRUE),
      expand = c(0, 0)
    ) +
    expand_limits(y = y_bracket + 0.06 * y_max) +
    labs(x = "Value (Relative to Remaining Service Area Mean)",
         y = "Density") +
    scale_color_manual(
      values = c(lost = "red", remaining = "blue"),
      labels = c(lost = "Wetland Loss (1985–2021)",
                 remaining = "Remaining Wetland (2021)"),
      name   = ""
    ) +
    theme_bw(base_size = 14) +
    theme(
      legend.position = "bottom",
      legend.title = element_text(size = 14),
      legend.text = element_text(size = 12),
      axis.title = element_text(size = 16),
      axis.text = element_text(size = 14),
      plot.title = element_text(size = 16),
      panel.grid.major = element_blank(),
      panel.grid.minor = element_blank()
    )
  
  message("All done! Returning plots and data.")
  list(
    pdf_plot      = pdf_plot,
    data          = results_ecdf_full,
    p_value       = p_val,
    sig_stars     = sig_stars
  )
}

test = readRDS(pth("CONUS", "Extractions and Summaries", "full_summary.rds"))
names = unique(test$sa_id)
# names = names[1:6]

datadir = "ECDF Functions"
res <- compare_remaining_vs_lost(datadir, bank_filter = names)

print(res$pdf_plot)

# ------------------------------------------
#$$$$$$$$$$$$$$$$$ FIGURE 4 $$$$$$$$$$$$$$$$$$$$


# ──────────────────────────────────────────────────────────────
#        1.  paths
# ──────────────────────────────────────────────────────────────
base_dir     <- datadir

pool_dir     <- file.path(base_dir, "Pooled HUC8 ECDFs")
bank_dir      <- file.path(base_dir, "Bank ECDF Functions")
loss_dir     <- file.path(base_dir, "SA ECDF Functions - Loss")

pool_lookup <- bank_pool_map %>%                 # pool_id  ↔  sa_id
  distinct(pool_id, sa_id = bank_name)

sa_multi <- pool_lookup$sa_id                    # ≥2 banks
sa_loss  <- sub("^ecdf_fn_(.*)\\.rds$", "\\1",
                list.files(loss_dir))            # every SA that has loss ECDF
sa_solo  <- setdiff(sa_loss, sa_multi)           # exactly 1 bank

# master list of “units” to iterate over
unit_df <- tibble(
  unit_name = c(sa_multi, sa_solo),
  type      = c(rep("multi", length(sa_multi)),
                rep("solo",  length(sa_solo)))
)

# ──────────────────────────────────────────────────────────────
# 3.  Helper functions  
# ──────────────────────────────────────────────────────────────
get_ecdf_data <- function(ecdf_obj) {
  dat <- environment(ecdf_obj)$x
  if (length(dat) > 1e6) dat <- sample(dat, 1e6)
  dat
}


as_num_vec <- function(v) {
  v <- suppressWarnings(as.numeric(v))
  v[is.finite(v)]
}

signed_area_diff <- function(ecdf2, ecdf1, grid_size = 1000) {
  x1 <- get_ecdf_data(ecdf1)
  x2 <- get_ecdf_data(ecdf2)
  grid <- seq(min(c(x1, x2)), max(c(x1, x2)), length.out = grid_size)
  delta_x <- diff(range(grid)) / (grid_size - 1)
  sum(ecdf1(grid) - ecdf2(grid)) * delta_x
}

quantile_diff <- function(x1, x2,
                          probs = c(.25, .5, .75, .95, .99)) {
  setNames(quantile(x1, probs) - quantile(x2, probs),
           paste0("qdiff_", probs * 100))
}

make_row <- function(unit_name, label, ecdf1, ecdf2) {
  if (!is.function(ecdf1) || !is.function(ecdf2))
    return(tibble(unit = unit_name, comparison = label,
                  mean_diff = NA, area_diff = NA,
                  qdiff_25 = NA, qdiff_50 = NA, qdiff_75 = NA,
                  qdiff_95 = NA, qdiff_99 = NA))
  
  x1 <- get_ecdf_data(ecdf1);  x2 <- get_ecdf_data(ecdf2)
  qd <- quantile_diff(x1, x2)
  tibble(unit = unit_name, comparison = label,
         mean_diff = mean(x1) - mean(x2),
         area_diff = signed_area_diff(ecdf1, ecdf2),
         !!!qd)
}

# ──────────────────────────────────────────────────────────────
# 4.  ECDF comparison worker
# ──────────────────────────────────────────────────────────────
process_unit <- function(unit_name, unit_type,
                         bank_dir, pool_dir, loss_dir) {
  # ---- BANK or pooled BANK ECDF ----
  if (unit_type == "multi") {
    pid       <- pool_lookup$pool_id[pool_lookup$sa_id == unit_name][1]
    bank_path <- file.path(pool_dir,
                           paste0("ecdf_pool_vs_sa_", pid, ".rds"))
    if (!file.exists(bank_path))
      return(list(status = "skip",
                  reason = "missing pooled file",
                  unit   = unit_name))
    
    bank_obj <- readRDS(bank_path)
    bank_ecdf_units  <- bank_obj$bank_ecdf$housing_units
    bank_ecdf_values <- bank_obj$bank_ecdf$housing_value
  } else {  # solo
    bank_path <- file.path(bank_dir,
                           paste0("bank_ecdf_fn_", unit_name, ".rds"))
    if (!file.exists(bank_path))
      return(list(status = "skip",
                  reason = "missing bank ECDF",
                  unit   = unit_name))
    
    bank_obj <- readRDS(bank_path)
    bank_ecdf_units  <- bank_obj$housing_units
    bank_ecdf_values <- bank_obj$housing_value
  }
  
  # ---- LOSS ECDF ----
  loss_path <- file.path(loss_dir,
                         paste0("ecdf_fn_", unit_name, ".rds"))
  if (!file.exists(loss_path))
    return(list(status = "skip",
                reason = "missing loss ECDF",
                unit   = unit_name))
  
  loss_obj <- readRDS(loss_path)
  loss_ecdf_units  <- loss_obj$housing_units
  loss_ecdf_values <- loss_obj$housing_value
  
  # ---- Raw vectors (numeric, finite) ----
  bank_units_vec <- as_num_vec(get_ecdf_data(bank_ecdf_units))
  bank_vals_vec  <- as_num_vec(get_ecdf_data(bank_ecdf_values))
  loss_units_vec <- as_num_vec(get_ecdf_data(loss_ecdf_units))
  loss_vals_vec  <- as_num_vec(get_ecdf_data(loss_ecdf_values))
  
  # ---- Safety checks ----
  if (length(bank_units_vec) == 0 || length(loss_units_vec) == 0)
    return(list(status = "skip",
                reason = "no numeric housing-unit data",
                unit   = unit_name))
  
  if (length(bank_vals_vec)  == 0 || length(loss_vals_vec)  == 0)
    return(list(status = "skip",
                reason = "no numeric housing-value data",
                unit   = unit_name))
  
  mean_units <- mean(bank_units_vec, na.rm = TRUE)
  mean_vals  <- mean(bank_vals_vec,  na.rm = TRUE)
  
  if (!is.finite(mean_units) || mean_units == 0 ||
      !is.finite(mean_vals)  || mean_vals  == 0)
    return(list(status = "skip",
                reason = "bank mean is NA or 0",
                unit   = unit_name))
  
  # ---- Normalize + down-sample ----
  loss_units_norm <- loss_units_vec / mean_units
  loss_vals_norm  <- loss_vals_vec  / mean_vals
  
  if (length(loss_units_norm) > 1e5)
    loss_units_norm <- sample(loss_units_norm, 1e5)
  if (length(loss_vals_norm)  > 1e5)
    loss_vals_norm  <- sample(loss_vals_norm,  1e5)
  
  ecdf_units_norm <- ecdf(loss_units_norm)
  ecdf_vals_norm  <- ecdf(loss_vals_norm)
  
  x_grid <- 10^seq(-2, 3.5, length.out = 750)
  
  tibble(
    unit       = unit_name,
    comparison = rep(c("normalized_loss_units",
                       "normalized_loss_values"), each = length(x_grid)),
    x          = rep(x_grid, 2),
    y          = c(ecdf_units_norm(x_grid),
                   ecdf_vals_norm(x_grid)),
    n_obs      = length(loss_units_norm)
  )
}

# ──────────────────────────────────────────────────────────────
# 5.  Run with progress bar
# ──────────────────────────────────────────────────────────────
plan(multisession, workers = 5)
handlers("txtprogressbar")

skip_log <- list()           # collects skipped units + reasons

with_progress({
  p <- progressor(along = nrow(unit_df))
  
  result_list <- future_lapply(seq_len(nrow(unit_df)), function(i) {
    p(unit_df$unit_name[i])
    out <- process_unit(unit_df$unit_name[i],
                        unit_df$type[i],
                        bank_dir, pool_dir, loss_dir)
    out
  })
})

# separate successful tibbles vs skip records
results_long_normalized <- bind_rows(
  result_list[vapply(result_list, is.data.frame, logical(1))]
)

fwrite(results_long_normalized,file=paste0(datadir,"\\poolednormalized_ecdfs_loss.csv"))

# Filter only housing unit rows
results_hu <- results_long_normalized %>%
  filter(comparison == "normalized_loss_units")

# Recalculate ecdf_summary for housing units only
ecdf_summary <- results_hu %>%
  group_by(x, comparison) %>%
  summarize(mean_y = weighted.mean(y, w = n_obs), .groups = "drop")
crossings <- ecdf_summary %>%
  # In case zero or negative x, remove them to avoid log-scale or approx issues
  filter(x > 0) %>%
  group_by(comparison) %>%
  summarize(
    # approx(..., rule=2) will clamp outside x-range rather than return NA
    below_1   = approx(x, mean_y, xout = 1,   rule = 2)$y,
    below_1_5 = approx(x, mean_y, xout = 1.5, rule = 2)$y,
    below_2   = approx(x, mean_y, xout = 2,   rule = 2)$y,
    .groups = "drop"
  ) %>%
  pivot_longer(
    cols = starts_with("below_"),
    names_to = "key",
    values_to = "cdf_value"
  ) %>%
  mutate(
    ref_x = case_when(
      key == "below_1"   ~ 1,
      key == "below_1_5" ~ 1.5,
      key == "below_2"   ~ 2
    ),
    percent_below = round(100 * cdf_value, 1),
    percent_above = round(100 * (1 - cdf_value), 1),
    label = paste0(percent_below, "% below\n", percent_above, "% above")
  )
# Assign colors to the reference lines
ref_line_colors <- c("gray10", "gray25", "gray40")
names(ref_line_colors) <- c("below_1", "below_1_5", "below_2")


# Merge colors into crossings
crossings <- crossings %>%
  mutate(line_color = ref_line_colors[key])
crossings_hu <- crossings %>% filter(comparison == "normalized_loss_units")
xmin_val <- min(results_long_normalized$x[results_long_normalized$x > 0], na.rm = TRUE)
fig3_ecdf <- ggplot(
  results_hu,
  aes(x = x, y = y, group = unit, lwd = n_obs)
) +
  geom_line(alpha = 0.2, aes(color = "Big Cypress Mitigation Bank")) +  # label the raw lines+
  geom_line(data = ecdf_summary_hu, aes(x = x, y = mean_y, color = "Mean ECDF"), inherit.aes = FALSE, linewidth = 1) +
  
  # Vertical dashed lines from xmin to crossing point only
  geom_segment(
    data = crossings_hu,
    aes(x = ref_x, xend = ref_x, y = 0, yend = cdf_value, color = key),
    linetype = "dashed",
    linewidth = 0.5,
    inherit.aes = FALSE
  ) +
  
  # Horizontal dashed lines from xmin to vertical line, color-coded
  geom_segment(
    data = crossings_hu,
    aes(x = xmin_val, xend = ref_x, y = cdf_value, yend = cdf_value, color = key),
    linetype = "dashed",
    linewidth = 0.5,
    inherit.aes = FALSE
  ) +
  # annotate("text",
  #          x = big_cypress_label_point$x + 50,
  #          y = big_cypress_label_point$y + 0.01,
  #          label = "Big Cypress\nMitigation Bank",
  #          color = "black",
  #          fontface = "bold",
  #          size = 4,
  #          hjust = 0.3) +
  # annotate("segment",
  #          x = big_cypress_label_point$x + 26,  # moved from +20 to +23
  #          xend = big_cypress_label_point$x + 6,  # was just big_cypress_label_point$x
  #          y = big_cypress_label_point$y + 0.01,
  #          yend = big_cypress_label_point$y,
  #          color = "black",
  #          arrow = arrow(length = unit(0.25, "cm")))+
  
  # Color mapping for reference lines
  scale_color_manual(
    values = c(
      "below_1" = "gray10",
      "below_1_5" = "gray25",
      "below_2" = "gray40",
      "Mean ECDF" = "black",
      "Big Cypress Mitigation Bank" = "#F8766D"
    ),
    guide = "none"  # remove legend
  )+
  
  labs(
    y = "Cumulative Density",
    x = "Value Relative to Bank Mean"
  ) +
  scale_x_log10(
    breaks = c(0.1, 1, 1.5, 2, 10, 100, 1000),
    labels = c("0.1", "1", "1.5", "2", "10", "100", "1000"),
    expand = c(0, 0)
  )+
  scale_y_continuous(expand = c(0, 0))+
  theme_bw(base_size = 14) +  # make all text larger
  theme(
    legend.position = "bottom",
    legend.title = element_text(size = 14),
    legend.text = element_text(size = 12),
    strip.text = element_blank(),  # remove facet title
    strip.background = element_blank(),
    axis.title = element_text(size = 16),
    axis.text = element_text(size = 14),
    panel.grid.major = element_blank(),
    panel.grid.minor = element_blank()
  ) +
  scale_linewidth_continuous(guide = "none")

fig3_ecdf
