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
  library(readr)
})


pth <- function(...) file.path("C:/Users/indumati/Box/Paper2_final", ...)

#-----------------------------------------------------#
# Load banks and service areas
#-----------------------------------------------------#
setwd("C:/Users/indumati/Box/Paper2_final")

serviceareas = readRDS(pth("Service Areas", "ServiceAreas_agg.rds"))
banks=readRDS(pth("Bank Footprints", "footprints_and_buffers.rds"))

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

#######


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
#  bbox_duplicates is already in memory; if not, readRDS it here
bank_pool_map <- bbox_duplicates %>% 
  mutate(pool_id = paste0("Pool_", row_number())) %>% 
  select(pool_id, identical_sa_ids) %>% 
  separate_rows(identical_sa_ids, sep = ",\\s*") %>% 
  rename(bank_name_raw = identical_sa_ids) %>% 
  mutate(bank_name = sanitize_name(bank_name_raw))

# ------------------------------------------------------------
#            2.  locate extract files on disk  (once)
# ------------------------------------------------------------

bank_extract_dir <- pth("Extractions and Summaries", "Extract Banks")
sa_extract_dir   <- pth("Extractions and Summaries", "Extract Service Areas")

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
                           paste0("ecdf_fn_", pid, ".rds")))
  invisible(res)   # keeps memory light
})

message("Done!  Results saved in: ", output_dir)

#----------------------------------
# %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

# Make PDF plot

# Separate the data processing from plotting for compare_bank_vs_lost
process_bank_vs_lost_data <- function(
    datadir,
    sample_size = 10000,
    seed = 123,
    bank_filter = NULL,
    bank_pool_map = NULL,
    pooled_dir = file.path(datadir, "Pooled HUC8 ECDFs"),
    loss_pick_strategy = c("largest_n", "first")
) {
  suppressPackageStartupMessages({
    library(tidyverse)
    library(scales)
    library(vctrs)
  })
  loss_pick_strategy <- match.arg(loss_pick_strategy)
  
  # Helpers ----
  read_ecdf_vals <- function(ecdf_obj) {
    if (is.null(ecdf_obj)) return(list(x = numeric(0), n = 0L))
    xs <- knots(ecdf_obj)
    n  <- tryCatch({ environment(ecdf_obj)$n }, error = function(e) length(xs))
    list(x = xs, n = n %||% length(xs))
  }
  
  pick_one_loss_for_pool <- function(pool_banks) {
    # choose one loss ECDF (prefer the one with largest N)
    candidates <- map(pool_banks, function(b) {
      loss_file <- file.path(datadir, "SA ECDF Functions - Loss", sprintf("ecdf_fn_%s.rds", b))
      if (!file.exists(loss_file)) return(NULL)
      ec <- readRDS(loss_file)$housing_units
      vals <- read_ecdf_vals(ec)
      if (length(vals$x) == 0) return(NULL)
      tibble(bank = b, file = loss_file, n = vals$n)
    }) %>% list_rbind()
    if (nrow(candidates) == 0) return(NULL)
    chosen <- if (loss_pick_strategy == "largest_n") {
      candidates |> arrange(desc(n)) |> slice(1)
    } else {
      candidates |> arrange(bank) |> slice(1)
    }
    chosen$file
  }
  
  message("Scanning directory for ECDF file names …")
  # These filenames drive which "bank names" exist
  sa_files <- list.files(file.path(datadir, "SA ECDF Functions"), full.names = FALSE)
  bank_names <- sub("^ecdf_fn_(.*)\\.rds$", "\\1", sa_files)
  
  if (!is.null(bank_filter)) {
    bank_names <- bank_names[bank_names %in% bank_filter]
  }
  
  # If pooling, split the names into pooled vs non-pooled
  pooled_present <- !is.null(bank_pool_map) && nrow(bank_pool_map) > 0
  pool_map <- NULL
  banks_in_pools <- character(0)
  pool_ids <- character(0)
  
  if (pooled_present) {
    # ensure we only consider pools that intersect the working bank_names
    pool_map <- bank_pool_map |>
      mutate(bank_name = bank_name %||% bank_name_raw) |>
      filter(bank_name %in% bank_names) |>
      select(pool_id, bank_name) |>
      distinct()
    banks_in_pools <- unique(pool_map$bank_name)
    pool_ids <- unique(pool_map$pool_id)
  }
  
  # Worklist
  solo_bank_names <- setdiff(bank_names, banks_in_pools)
  total_tasks <- length(solo_bank_names) + length(pool_ids)
  if (total_tasks == 0) stop("No ECDF files found for specified banks in directory: ", datadir)
  
  results_list <- vector("list", total_tasks * 2)  # each task adds 'remaining' and 'lost'
  j <- 1
  successful_tasks <- 0
  set.seed(seed)
  
  # ---------- Process non-pooled banks (original behavior) ----------
  for (i in seq_along(solo_bank_names)) {
    bank <- solo_bank_names[i]
    message(sprintf("[Solo %d/%d] Processing bank: %s", i, total_tasks, bank))
    
    sa_file <- file.path(datadir, "SA ECDF Functions",        sprintf("ecdf_fn_%s.rds", bank))
    loss_file <- file.path(datadir, "SA ECDF Functions - Loss",   sprintf("ecdf_fn_%s.rds",      bank))
    
    if (!file.exists(sa_file) || !file.exists(loss_file)) {
      message("  - Skipping – missing ECDF file(s)")
      next
    }
    
    bank_ecdf  <- readRDS(sa_file)$housing_units
    loss_ecdf  <- readRDS(loss_file)$housing_units
    if (is.null(bank_ecdf) || is.null(loss_ecdf)) {
      message("  - Skipping – housing_units ECDF missing")
      next
    }
    
    bvals <- read_ecdf_vals(bank_ecdf)
    lvals <- read_ecdf_vals(loss_ecdf)
    if (length(bvals$x) == 0 || length(lvals$x) == 0) {
      message("  - Skipping – ECDF contains no data")
      next
    }
    
    bank_mean <- mean(bvals$x)
    
    # Downsample (values only; keep original N for weights)
    bxs <- if (length(bvals$x) > sample_size) {
      message(sprintf("  - Downsampled remaining: %s -> %s", comma(bvals$n), comma(sample_size)))
      sample(bvals$x, sample_size)
    } else bvals$x
    
    lxs <- if (length(lvals$x) > sample_size) {
      message(sprintf("  - Downsampled lost: %s -> %s", comma(lvals$n), comma(sample_size)))
      sample(lvals$x, sample_size)
    } else lvals$x
    
    row_weight <- lvals$n / length(lxs)
    
    remaining <- tibble(
      bank = bank, status = "bank", value = bxs,
      normalized_value = bxs / bank_mean,
      weight = row_weight
    ) |> filter(is.finite(normalized_value) & normalized_value > 0)
    
    lost <- tibble(
      bank = bank, status = "lost", value = lxs,
      normalized_value = lxs / bank_mean,
      weight = row_weight
    ) |> filter(is.finite(normalized_value) & normalized_value > 0)
    
    results_list[[j]] <- remaining; j <- j + 1
    results_list[[j]] <- lost;      j <- j + 1
    successful_tasks  <- successful_tasks + 1
  }
  
  # ---------- Process pooled service areas ----------
  if (pooled_present && length(pool_ids) > 0) {
    for (k in seq_along(pool_ids)) {
      pid <- pool_ids[k]
      message(sprintf("[Pool %d/%d] Processing pool: %s", length(solo_bank_names) + k, total_tasks, pid))
      
      pooled_file <- file.path(pooled_dir, sprintf("ecdf_pool_vs_sa_%s.rds", pid))
      if (!file.exists(pooled_file)) {
        message("  - Skipping – pooled ECDF file not found: ", pooled_file)
        next
      }
      pooled_obj <- readRDS(pooled_file)
      
      # Remaining service area ECDF (pooled)
      sa_ecdf <- pooled_obj$sa_rem_ecdf$housing_units
      bank_ecdf_pooled <- pooled_obj$bank_ecdf$housing_units  # for normalization baseline
      
      if (is.null(sa_ecdf) || is.null(bank_ecdf_pooled)) {
        message("  - Skipping – pooled housing_units ECDF missing")
        next
      }
      
      sa_vals   <- read_ecdf_vals(sa_ecdf)
      bank_vals <- read_ecdf_vals(bank_ecdf_pooled)
      if (length(sa_vals$x) == 0 || length(bank_vals$x) == 0) {
        message("  - Skipping – pooled ECDF has no data")
        next
      }
      
      bank_mean <- mean(bank_vals$x)
      
      # Downsample remaining (pooled SA)
      sa_xs <- if (length(sa_vals$x) > sample_size) {
        message(sprintf("  - Downsampled pooled remaining: %s -> %s", comma(sa_vals$n), comma(sample_size)))
        sample(sa_vals$x, sample_size)
      } else sa_vals$x
      
      # Choose ONE loss ECDF from the pool's member banks
      pool_banks <- pool_map |> filter(pool_id == pid) |> pull(bank_name) |> unique()
      chosen_loss_file <- pick_one_loss_for_pool(pool_banks)
      if (is.null(chosen_loss_file)) {
        message("  - Skipping – no usable loss ECDF in this pool")
        next
      }
      loss_ecdf <- readRDS(chosen_loss_file)$housing_units
      lvals <- read_ecdf_vals(loss_ecdf)
      if (length(lvals$x) == 0) {
        message("  - Skipping – chosen loss ECDF has no data")
        next
      }
      lxs <- if (length(lvals$x) > sample_size) {
        message(sprintf("  - Downsampled pooled lost: %s -> %s", comma(lvals$n), comma(sample_size)))
        sample(lvals$x, sample_size)
      } else lvals$x
      
      row_weight <- lvals$n / length(lxs)
      
      # Use pool_id as the series label
      remaining <- tibble(
        bank = pid, status = "bank", value = sa_xs,
        normalized_value = sa_xs / bank_mean,
        weight = row_weight
      ) |> filter(is.finite(normalized_value) & normalized_value > 0)
      
      lost <- tibble(
        bank = pid, status = "lost", value = lxs,
        normalized_value = lxs / bank_mean,
        weight = row_weight
      ) |> filter(is.finite(normalized_value) & normalized_value > 0)
      
      results_list[[j]] <- remaining; j <- j + 1
      results_list[[j]] <- lost;      j <- j + 1
      successful_tasks  <- successful_tasks + 1
    }
  }
  
  if (successful_tasks == 0) stop("No banks/pools with both remaining and loss ECDFs were processed.")
  message("Successfully processed ", successful_tasks, " series (banks or pools).")
  
  results <- bind_rows(results_list)
  rm(results_list); gc()
  
  message("Computing ECDF data …")
  results_ecdf_full <- results %>%
    group_by(bank, status) %>%
    arrange(normalized_value) %>%
    mutate(ecdf_y = row_number() / n()) %>%
    ungroup() %>%
    mutate(log_value = log10(normalized_value))
  
  # Return just the processed data
  return(results_ecdf_full)
}

# Separate plotting function that takes processed data
create_bank_vs_lost_plot <- function(data, 
                                     plot_title = NULL,
                                     x_limits = c(0.01, 100),
                                     x_label = "Flood Protection Index (Relative to Remaining Service Area Mean)",
                                     y_label = "Density",
                                     colors = c(lost = "red", bank = "blue"),
                                     show_stats = TRUE,
                                     base_size = 14) {
  
  # Calculate statistics
  loss_mean <- weighted.mean(
    data$normalized_value[data$status == "lost"],
    w = data$weight[data$status == "lost"], na.rm = TRUE
  )
  
  # t-test (unweighted)
  p_val <- t.test(normalized_value ~ status, data = data)$p.value
  sig_stars <- case_when(
    p_val < 0.001 ~ "***",
    p_val < 0.01  ~ "**",
    p_val < 0.05  ~ "*",
    TRUE          ~ "ns"
  )
  
  # Calculate bracket position - use unweighted density for consistent scaling
  x_bank <- data$normalized_value[data$status == "bank"]
  x_lost <- data$normalized_value[data$status == "lost"]
  
  dens_rem  <- density(x_bank)
  dens_lost <- density(x_lost)
  y_max     <- max(dens_rem$y, dens_lost$y, na.rm = TRUE)
  y_bracket <- y_max * .3
  y_tick    <- y_max * 0.3
  x_left  <- min(1, loss_mean)
  x_right <- max(1, loss_mean)
  
  # Create base plot
  p <- ggplot(data, aes(x = normalized_value, colour = status, weight = weight)) +
    stat_density(geom = "line", position = "identity", linewidth = 1.2) +
    geom_vline(xintercept = 1,         linetype = "dashed", colour = "blue", linewidth = .8) +
    geom_vline(xintercept = loss_mean, linetype = "dashed", colour = "red",  linewidth = .8) +
    scale_x_log10(
      limits = x_limits,
      breaks = sort(unique(c(0.1, 1, 10, 100, signif(loss_mean, 2)))),
      labels = function(x) format(x, digits = 2, trim = TRUE),
      expand = c(0, 0)
    ) +
    labs(x = x_label,
         y = y_label,
         title = plot_title) +
    scale_color_manual(
      values = colors,
      labels = c(lost = "Wetland Loss (1985–2021)",
                 bank = "Remaining Wetland (2021)"),
      name   = ""
    ) +
    theme_bw(base_size = base_size) +
    theme(
      legend.position = "bottom",
      legend.title = element_text(size = base_size),
      legend.text = element_text(size = base_size - 2),
      axis.title = element_text(size = base_size + 2),
      axis.text = element_text(size = base_size),
      plot.title = element_text(size = base_size + 2),
      panel.grid.major = element_blank(),
      panel.grid.minor = element_blank()
    )
  
  # Add statistical annotations if requested
  if (show_stats) {
    p <- p +
      annotate("segment", x = x_left,  xend = x_right, y = y_bracket, yend = y_bracket) +
      annotate("segment", x = x_left,  xend = x_left,  y = y_bracket, yend = y_tick) +
      annotate("segment", x = x_right, xend = x_right, y = y_bracket, yend = y_tick) +
      annotate("text", x = (x_left + x_right)/2 - 0.2, y = y_bracket + 0.03 * y_max,
               label = sig_stars, vjust = 0) +
      expand_limits(y = y_bracket + 0.06 * y_max)
  }
  
  return(list(
    plot = p,
    p_value = p_val,
    sig_stars = sig_stars,
    loss_mean = loss_mean
  ))
}


processed_data <- process_bank_vs_lost_data(
  datadir,
  bank_filter = names,
  bank_pool_map = bank_pool_map,
  pooled_dir = file.path(datadir, "Pooled HUC8 ECDFs"),
  loss_pick_strategy = "largest_n"
)

# Create plot
plot_result <- create_bank_vs_lost_plot(processed_data)
print(plot_result$plot)

ggsave(file.path("CONUS/Graphics/pdfplot_pooled.png"),  plot_result$plot, width=12, height=8, dpi=300)

# ------------------------------------------
#$$$$$$$$$$$$$$$$$ FIGURE 4 ECDF $$$$$$$$$$$$$$$$$$$$


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

# Pool representative workflow 
library(dplyr)
library(data.table)
library(future)
library(future.apply)
library(progressr)

# Setup directories and file lists
base_dir <- datadir
bankfiles <- list.files(file.path(base_dir, "Bank ECDF Functions"), full.names = TRUE)
lossfiles <- list.files(file.path(base_dir, "SA ECDF Functions - Loss"), full.names = TRUE)

# Get available unit names from bank files
bank_names <- list.files(file.path(base_dir, "Bank ECDF Functions"))
bank_names <- sub("^bank_ecdf_fn_(.*)/.rds$", "/1", bank_names)

loss_names <- list.files(file.path(base_dir, "SA ECDF Functions - Loss"))
loss_names <- sub("^ecdf_fn_(.*)/.rds$", "/1", loss_names)

# Units that have both bank AND loss files
complete_units <- intersect(bank_names, loss_names)
cat("Units with both bank and loss files:", length(complete_units), "\n")

# Create pool lookup
pool_lookup <- bank_pool_map %>%
  distinct(pool_id, sa_id = bank_name)

# Get ONE representative per pool from units that have complete files
pools_with_complete <- pool_lookup %>%
  filter(sa_id %in% complete_units) %>%
  group_by(pool_id) %>%
  slice(1) %>%  # Take first complete unit per pool
  ungroup()

# Get representative units for multi-bank pools
sa_multi <- pools_with_complete$sa_id

# Solo units are those NOT in any pool AND have complete files
sa_solo <- setdiff(complete_units, pool_lookup$sa_id)

# Final list of units to process
names <- c(sa_multi, sa_solo)

cat("Pool representatives:", length(sa_multi), "\n")
cat("Solo units:", length(sa_solo), "\n")
cat("Total units to process:", length(names), "\n")

# Verify coverage
pools_covered <- pools_with_complete$pool_id
total_pools <- unique(pool_lookup$pool_id)
cat("Pools covered:", length(pools_covered), "out of", length(total_pools), "\n")

# Setup x_grid (same as your original)
x_grid <- 10^seq(-2, 3.5, length.out = 750)

process_unit <- function(unit_name, datadir=datadir, bankfiles=bankfiles, lossfiles=lossfiles, x_grid=x_grid) {
  # Always load loss data
  lossfile <- paste0(datadir, "/SA ECDF Functions - Loss/ecdf_fn_", unit_name, ".rds")
  loss_data <- if (lossfile %in% lossfiles) readRDS(lossfile) else NA
  
  # Determine if unit is in a pool
  pool_id <- pool_lookup$pool_id[pool_lookup$sa_id == unit_name]
  
  if (length(pool_id) > 0) {
    # POOLED UNIT: Load from pooled directory
    poolfile <- paste0(datadir, "/Pooled HUC8 ECDFs/ecdf_pool_vs_sa_", pool_id[1], ".rds")
    if (file.exists(poolfile)) {
      pool_data <- readRDS(poolfile)
      bank_data <- pool_data$bank_ecdf
    } else {
      return(NULL)
    }
  } else {
    # SOLO UNIT: Load from individual bank directory
    bankfile <- paste0(datadir, "/Bank ECDF Functions/bank_ecdf_fn_", unit_name, ".rds")
    bank_data <- if (bankfile %in% bankfiles) readRDS(bankfile) else NA
  }
  
  # Validate required objects
  if (!is.list(bank_data) || !is.list(loss_data)) return(NULL)
  if (!is.function(bank_data$housing_units) || !is.function(loss_data$housing_units)) return(NULL)
  if (!is.function(bank_data$housing_value) || !is.function(loss_data$housing_value)) return(NULL)
  
  # Extract raw values
  bank_vals  <- environment(bank_data$housing_units)$x
  bank_hval  <- environment(bank_data$housing_value)$x
  loss_vals  <- environment(loss_data$housing_units)$x
  loss_hval  <- environment(loss_data$housing_value)$x
  
  # Compute means
  bank_mean_vals <- mean(bank_vals, na.rm = TRUE)
  bank_mean_hval <- mean(bank_hval, na.rm = TRUE)
  
  # Normalize accordingly
  loss_vals_norm <- loss_vals / bank_mean_vals
  loss_hval_norm <- loss_hval / bank_mean_hval
  
  n_obs <- length(loss_vals_norm)
  
  # Optional downsampling
  if (n_obs > 1e5) {
    sampled_indices <- sample(seq_len(n_obs), 1e5)
    loss_vals_norm <- loss_vals_norm[sampled_indices]
    loss_hval_norm <- loss_hval_norm[sampled_indices]
  }
  
  # Evaluate ECDFs
  ecdf_vals  <- ecdf(loss_vals_norm)
  ecdf_hvals <- ecdf(loss_hval_norm)
  
  # Output both sets as tibbles
  result_units <- tibble(
    unit = unit_name,
    comparison = "normalized_loss_units",
    x = x_grid,
    y = ecdf_vals(x_grid),
    n_obs = n_obs
  )
  
  result_values <- tibble(
    unit = unit_name,
    comparison = "normalized_loss_values",
    x = x_grid,
    y = ecdf_hvals(x_grid),
    n_obs = n_obs
  )
  
  # Return combined result
  bind_rows(result_units, result_values)
}

# Your exact parallel processing workflow
plan(multisession, workers= 10)
handlers("txtprogressbar")

with_progress({
  p <- progressor(along = names)
  
  results_long <- future_lapply(names, function(n) {
    p(message = n)
    process_unit(n, datadir, bankfiles, lossfiles, x_grid)
  })
})

results_long <- Filter(Negate(is.null), results_long)
results_long_normalized <- bind_rows(results_long)


fwrite(results_long_normalized,file=paste0(datadir,"\\poolednormalized_ecdfs_loss.csv"))




# Highest number of loss cells:
results_long_normalized = read.csv("ECDF Functions\\poolednormalized_ecdfs_loss.csv")
# Filter only housing unit rows
results_hu <- results_long_normalized %>%
  filter(comparison == "normalized_loss_units")

# Recalculate ecdf_summary for housing units only
ecdf_summary <- results_hu %>%
  group_by(x, comparison) %>%
  summarize(mean_y = weighted.mean(y, w = n_obs), .groups = "drop")
crossings <- ecdf_summary %>%
  # In case you have zero or negative x, remove them to avoid log-scale or approx issues
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
    # Example label: "30% below\n70% above"
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
f_mean <- approxfun(ecdf_summary_hu$x, ecdf_summary_hu$mean_y, rule = 2, ties = mean)

# Reference x’s (clamped to the curve’s domain for safety on log scale)
ref_tbl <- tibble(
  key   = c("below_1", "below_1_5", "below_2"),
  ref_x = c(1, 1.5, 2)
) %>%
  mutate(
    ref_x = pmax(min(ecdf_summary_hu$x, na.rm = TRUE),
                 pmin(ref_x, max(ecdf_summary_hu$x, na.rm = TRUE)))
  )
crossings_hu <- ref_tbl %>%
  mutate(
    cdf_value     = f_mean(ref_x),
    percent_below = round(100 * cdf_value, 1),
    percent_above = round(100 * (1 - cdf_value), 1),
    label         = paste0(percent_below, "% below\n", percent_above, "% above"),
    line_color    = c("below_1" = "gray10", "below_1_5" = "gray25", "below_2" = "gray40")[key]
  )
xmin_val <- min(ecdf_summary_hu$x, na.rm = TRUE)
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

ggsave(file.path("CONUS/Graphics/ecdfplot_pooled.png"),  fig3_ecdf, width=12, height=8, dpi=300)


################# Get actual numbers

ecdf_csv <- file.path(datadir, "poolednormalized_ecdfs_loss.csv")

results_long_normalized <- read_csv(ecdf_csv, show_col_types = FALSE)

# Keep housing units ECDFs that are normalized by bank mean
results_hu <- results_long_normalized %>%
  filter(comparison == "normalized_loss_units")

# Observation-weighted mean ECDF (y is already a CDF value)
ecdf_summary_hu <- results_hu %>%
  group_by(x) %>%
  summarize(mean_y = weighted.mean(y, w = n_obs, na.rm = TRUE), .groups = "drop")

# Turn into a function (clamps outside the domain)
f_mean <- approxfun(ecdf_summary_hu$x, ecdf_summary_hu$mean_y, rule = 2, ties = mean)

# Evaluate the percent ABOVE thresholds
cuts <- c(`>1x` = 1, `>1.5x` = 1.5, `>2x` = 2, `>10x` = 10)
percent_above_tbl <- tibble(
  cut = names(cuts),
  x   = as.numeric(cuts),
  cdf = f_mean(x),
  pct_above = (1 - cdf) * 100
) %>%
  mutate(pct_above = round(pct_above, 1))

percent_above_tbl

cut_pts <- c(`gt1` = 1, `gte1.5` = 1.5, `gte2` = 2, `gte4.5` = 4.5,  `gt10` = 10)

# helper function: 1 - ECDF(threshold)
tail_share <- function(df, threshold) {
  1 - approx(x = df$x, y = df$y, xout = threshold, rule = 2)$y
}

# compute tail shares per unit
shares_per_unit <- results_long_normalized %>%
  filter(comparison == "normalized_loss_units") %>%
  group_by(unit) %>%
  group_modify(~{
    unit_n <- max(.x$n_obs, na.rm = TRUE)
    shares <- map_dbl(cut_pts, function(thresh) tail_share(.x, thresh))
    tibble(n_obs = unit_n, !!!setNames(shares * 100, names(cut_pts)))
  }) %>%
  ungroup()
# overall weighted mean
overall_wtd <- shares_per_unit %>%
  select(-unit) %>%
  summarise(across(-n_obs, ~ weighted.mean(.x, w = n_obs, na.rm = TRUE)))
overall_wtd


datadir          <- "L:/Wetland Flood Mitigation/ECDF Functions"
bank_ecdf_dir    <- file.path(datadir, "Bank ECDF Functions")
sa_loss_ecdf_dir <- file.path(datadir, "SA ECDF Functions - Loss")
pooled_dir       <- file.path(datadir, "Pooled HUC8 ECDFs")

# bank_pool_map already in memory 

# Pool lookup: pool_id per service-area id 
pool_lookup <- bank_pool_map %>% distinct(pool_id, sa_id = bank_name)

# Build unit list (one rep per pool + all solo) 
if (!exists("names")) {
  bank_units  <- sub("^bank_ecdf_fn_(.*)\\.rds$", "\\1", list.files(bank_ecdf_dir))
  loss_units  <- sub("^ecdf_fn_(.*)\\.rds$", "\\1", list.files(sa_loss_ecdf_dir))
  complete    <- intersect(bank_units, loss_units)
  pools_with_complete <- pool_lookup %>%
    filter(sa_id %in% complete) %>%
    group_by(pool_id) %>% slice(1) %>% ungroup()
  sa_multi <- pools_with_complete$sa_id
  sa_solo  <- setdiff(complete, pool_lookup$sa_id)
  names    <- c(sa_multi, sa_solo)
}

# Helper to get the correct BANK ecdf object (pooled vs solo)
get_bank_ecdf <- function(unit_name) {
  pool_id <- pool_lookup$pool_id[pool_lookup$sa_id == unit_name]
  if (length(pool_id) > 0) {
    poolfile <- file.path(pooled_dir, paste0("ecdf_pool_vs_sa_", pool_id[1], ".rds"))
    if (!file.exists(poolfile)) return(NULL)
    readRDS(poolfile)$bank_ecdf
  } else {
    bankfile <- file.path(bank_ecdf_dir, paste0("bank_ecdf_fn_", unit_name, ".rds"))
    if (!file.exists(bankfile)) return(NULL)
    readRDS(bankfile)
  }
}

# Core computation: mean(loss) / mean(bank) for the HOUSING-UNITS index
ratio_df <- map_dfr(names, function(nm) {
  loss_pth <- file.path(sa_loss_ecdf_dir, paste0("ecdf_fn_", nm, ".rds"))
  if (!file.exists(loss_pth)) return(NULL)
  
  bank_ecdf <- get_bank_ecdf(nm)
  if (is.null(bank_ecdf)) return(NULL)
  
  # Pull raw x vectors from ECDF closures (housing_units index)
  if (!is.function(bank_ecdf$housing_units)) return(NULL)
  loss_ecdf <- readRDS(loss_pth)
  if (!is.function(loss_ecdf$housing_units)) return(NULL)
  
  bank_vals <- environment(bank_ecdf$housing_units)$x
  loss_vals <- environment(loss_ecdf$housing_units)$x
  
  # Guardrails
  if (!is.numeric(bank_vals) || !is.numeric(loss_vals)) return(NULL)
  if (!length(bank_vals) || !length(loss_vals))       return(NULL)
  
  tibble(
    unit      = nm,
    mean_bank = mean(bank_vals, na.rm = TRUE),
    mean_loss = mean(loss_vals, na.rm = TRUE),
    n_pixels  = length(loss_vals)
  )
}) %>%
  filter(is.finite(mean_bank), mean_bank > 0,
         is.finite(mean_loss), mean_loss > 0) %>%
  mutate(ratio_mean = mean_loss / mean_bank)

# Area-weighted “average across all service areas”
area_weighted_avg <- with(ratio_df, weighted.mean(ratio_mean, w = n_pixels, na.rm = TRUE))
max_row <- ratio_df %>% arrange(desc(ratio_mean)) %>% slice(1)

# Optional extras if you want them for robustness notes
p50_ratio <- median(ratio_df$ratio_mean, na.rm = TRUE)
p95_ratio <- quantile(ratio_df$ratio_mean, 0.95, na.rm = TRUE)

cat(glue(
  "wetlands lost to development exist near prior developed areas and therefore ",
  "provide relatively high levels of downstream flood protection compared to ",
  "wetlands created in compensation--on average \u2248 {round(area_weighted_avg, 1)} times as much, ",
  "though in some cases up to {round(max_row$ratio_mean, 0)} times."
))

max_row
