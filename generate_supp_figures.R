datadir <- Sys.getenv("DATA_DIR", "C:/Users/indumati/Box/Paper2_final")

##########################
######## FIGURE S1 #######
##########################

# ECDFs -----

library(ggplot2)
library(tidyr)
library(patchwork)
library(purrr)
library(scales)

# =========================
# 1) Update process_unit to also handle housing_value & population
# =========================
process_unit <- function(unit_name, datadir=datadir, bankfiles=bankfiles, lossfiles=lossfiles, x_grid=x_grid) {
  # Always load loss data
  lossfile <- file.path(datadir, "SA ECDF Functions - Loss", paste0("ecdf_fn_", unit_name, ".rds"))
  loss_data <- if (lossfile %in% lossfiles) readRDS(lossfile) else NA
  
  # Determine if unit is in a pool
  pool_id <- pool_lookup$pool_id[pool_lookup$sa_id == unit_name]
  
  if (length(pool_id) > 0) {
    poolfile <- file.path(datadir, "Pooled HUC8 ECDFs", paste0("ecdf_pool_vs_sa_", pool_id[1], ".rds"))
    if (file.exists(poolfile)) {
      pool_data <- readRDS(poolfile)
      bank_data <- pool_data$bank_ecdf
    } else {
      return(NULL)
    }
  } else {
    # SOLO UNIT: Load  individual bank directory
    bankfile <- file.path(datadir, "Bank ECDF Functions", paste0("bank_ecdf_fn_", unit_name, ".rds"))
    bank_data <- if (bankfile %in% bankfiles) readRDS(bankfile) else NA
  }
  
  need_lists <- is.list(bank_data) && is.list(loss_data)
  need_funcs <- all(
    is.function(bank_data$housing_units),
    is.function(loss_data$housing_units),
    is.function(bank_data$housing_value),
    is.function(loss_data$housing_value)
  )
  if (!need_lists || !need_funcs) return(NULL)
  
  # Extract raw values
  bank_units <- environment(bank_data$housing_units)$x
  bank_value <- environment(bank_data$housing_value)$x
  loss_units <- environment(loss_data$housing_units)$x
  loss_value <- environment(loss_data$housing_value)$x
  
  have_pop <- is.function(bank_data$population) && is.function(loss_data$population)
  if (have_pop) {
    bank_pop <- environment(bank_data$population)$x
    loss_pop <- environment(loss_data$population)$x
  }
  
  # Compute means
  bank_mean_units <- mean(bank_units, na.rm = TRUE)
  bank_mean_value <- mean(bank_value, na.rm = TRUE)
  bank_mean_pop   <- if (have_pop) mean(bank_pop, na.rm = TRUE) else NA_real_
  
  # Normalize accordingly
  loss_units_norm <- loss_units / bank_mean_units
  loss_value_norm <- loss_value / bank_mean_value
  loss_pop_norm   <- if (have_pop) loss_pop / bank_mean_pop else NULL
  
  # Downsamplign
  n_obs <- length(loss_units_norm)
  if (n_obs > 1e5) {
    sampled_indices <- sample(seq_len(n_obs), 1e5)
    loss_units_norm <- loss_units_norm[sampled_indices]
    loss_value_norm <- loss_value_norm[sampled_indices]
    if (have_pop) loss_pop_norm <- loss_pop_norm[sampled_indices]
    n_obs <- length(sampled_indices)
  }
  
  # Evaluate ECDFs
  ecdf_units <- ecdf(loss_units_norm)
  ecdf_value <- ecdf(loss_value_norm)
  
  out_units <- tibble(
    unit = unit_name,
    comparison = "normalized_loss_units",
    x = x_grid,
    y = ecdf_units(x_grid),
    n_obs = n_obs
  )
  
  out_value <- tibble(
    unit = unit_name,
    comparison = "normalized_loss_values",
    x = x_grid,
    y = ecdf_value(x_grid),
    n_obs = n_obs
  )
  
  if (have_pop) {
    ecdf_pop <- ecdf(loss_pop_norm)
    out_pop <- tibble(
      unit = unit_name,
      comparison = "normalized_loss_population",
      x = x_grid,
      y = ecdf_pop(x_grid),
      n_obs = n_obs
    )
    return(bind_rows(out_units, out_value, out_pop))
  } else {
    return(bind_rows(out_units, out_value))
  }
}

# =========================
# 2) Run parallel 
# =========================
plan(multisession, workers = 10)
handlers("txtprogressbar")

with_progress({
  p <- progressor(along = names)
  results_long <- future_lapply(names, function(nm) {
    p(message = nm)
    process_unit(nm, datadir, bankfiles, lossfiles, x_grid)
  })
})

results_long <- Filter(Negate(is.null), results_long)
results_long_normalized <- bind_rows(results_long)

# outfile <- file.path(datadir, "normalized_ecdfs_pophv.csv")
# fwrite(results_long_normalized, file = outfile)

#####################################


make_ecdf_plot <- function(df_all, comparison_key, panel_title = NULL,
                           ref_xs = c(1, 2),
                           raw_color = "#F8766D",        
                           mean_color = "black") {        
  # Filter to the metric
  df <- df_all %>% filter(comparison == comparison_key)
  
  # Weighted mean ECDF across units
  ecdf_summary <- df %>%
    group_by(x, comparison) %>%
    summarize(mean_y = weighted.mean(y, w = n_obs), .groups = "drop")
  
  # Safe approx function and x-domain clamps
  f_mean <- approxfun(ecdf_summary$x, ecdf_summary$mean_y, rule = 2, ties = mean)
  xmin_val <- min(ecdf_summary$x, na.rm = TRUE)
  xmax_val <- max(ecdf_summary$x, na.rm = TRUE)
  
  ref_tbl <- tibble(
    key   = paste0("below_", ref_xs),
    ref_x = pmax(xmin_val, pmin(ref_xs, xmax_val))
  ) %>%
    mutate(
      cdf_value     = f_mean(ref_x),
      percent_below = round(100 * cdf_value, 1),
      percent_above = round(100 * (1 - cdf_value), 1),
      label         = paste0(percent_below, "% below\n", percent_above, "% above")
    )
  
  # Line colors 
  ref_colors <- setNames(c("gray10", "gray25", "gray40"), paste0("below_", ref_xs))
  color_map <- c(
    ref_colors,
    "Mean ECDF" = mean_color,
    "Big Cypress Mitigation Bank" = raw_color
  )
  
  ggplot(
    df,
    aes(x = x, y = y, group = unit, lwd = n_obs)
  ) +
    geom_line(alpha = 0.2, aes(color = "Big Cypress Mitigation Bank")) +
    geom_line(
      data = ecdf_summary,
      aes(x = x, y = mean_y, color = "Mean ECDF"),
      inherit.aes = FALSE,
      linewidth = 1
    ) +
    geom_segment(
      data = ref_tbl,
      aes(x = ref_x, xend = ref_x, y = 0, yend = cdf_value, color = key),
      linetype = "dashed",
      linewidth = 0.5,
      inherit.aes = FALSE
    ) +
    geom_segment(
      data = ref_tbl,
      aes(x = xmin_val, xend = ref_x, y = cdf_value, yend = cdf_value, color = key),
      linetype = "dashed",
      linewidth = 0.5,
      inherit.aes = FALSE
    ) +
    scale_color_manual(values = color_map, guide = "none") +
    labs(
      title = panel_title %||% comparison_key,
      y = "Cumulative Density",
      x = "Flood Protection Index of Developed Wetlands, Relative to Bank Mean"
    ) +
    scale_x_log10(
      breaks = c(0.1, 1, 2, 10, 100, 1000),
      labels = c("0.1", "1", "2", "10", "100", "1000"),
      expand = c(0, 0)
    ) +
    scale_y_continuous(expand = c(0, 0)) +
    theme_bw(base_size = 14) +
    theme(
      plot.title = element_text(hjust =0.5),
      legend.position = "bottom",
      legend.title = element_text(size = 14),
      legend.text = element_text(size = 12),
      strip.text = element_blank(),
      strip.background = element_blank(),
      axis.title = element_text(size = 16),
      axis.text = element_text(size = 14),
      panel.grid.major = element_blank(),
      panel.grid.minor = element_blank()
    ) +
    scale_linewidth_continuous(guide = "none")
}

# =========================
# 4) Build the two panels with different palettes and display side-by-side
# =========================
p_value <- make_ecdf_plot(
  results_long_normalized,
  comparison_key = "normalized_loss_values",
  panel_title = "Housing Value",
  raw_color = "#9ecae1",   
  mean_color = "#1f78b4"   
)

p_pop <- make_ecdf_plot(
  results_long_normalized,
  comparison_key = "normalized_loss_population",
  panel_title = "Population",
  raw_color = "#b2df8a",   
  mean_color = "#33a02c"   
)+
  labs(y = NULL) +                                 # remove y-axis title
  theme(axis.title.y = element_blank())           

fig3_value_pop <- p_value | p_pop

##   remove x-axis titles inside each panel
p_value_no_x <- p_value + labs(x = NULL)
p_pop_no_x   <- p_pop   + labs(x = NULL)

##   combine the two panels
fig3_value_pop <- p_value_no_x | p_pop_no_x

##   add ONE global x-axis label via plot_annotation()
fig3_value_pop <- fig3_value_pop +
  plot_annotation(
    caption = "Flood Protection Index of Developed Wetlands, Relative to Bank Mean",
    theme = theme(
      plot.caption      = element_text(size = 16, hjust = 0.5, vjust = -0.5),
      axis.title.x      = element_blank(),   # ensure inner titles stay blank
      plot.margin       = margin(b = 20)     # a little space for the caption
    )
  )

fig3_value_pop
# Add shared x-axis label
fig3_value_pop <- add_sub(fig3_value_pop,
                          "Flood Protection Index of Developed Wetlands, Relative to Bank Mean",
                          vjust = 1, size = 12
)
fig3_value_pop

####################################
############ PDF ###################
####################################

`%||%` <- rlang::`%||%`   

# ──────────────────────────────────────────────────────────────────────────────
# 1)  General-purpose comparison function 
# ──────────────────────────────────────────────────────────────────────────────
compare_remaining_vs_lost_pooled <- function(datadir,
                                             bank_pool_map,
                                             bank_summary_path,
                                             metric = c("housing_units",
                                                        "housing_value",
                                                        "population"),
                                             sample_size = 10000,
                                             seed = 123,
                                             filter = NULL) {
  metric <- match.arg(metric)   # choose one of the three
  
  # ── load bank_summary & set up ────────────────────────────────────────────
  bank_summary <- readRDS(bank_summary_path)
  valid_sa_names <- unique(bank_summary$bank_name)   # column with SA names
  
  pool_representatives <- bank_pool_map %>%          # 1 rep per pool
    group_by(pool_id) %>% slice(1) %>% ungroup()
  
  pooled_banks <- unique(bank_pool_map$bank_name)
  
  sa_files <- list.files(file.path(datadir, "SA ECDF Functions"), full.names = FALSE)
  sa_names <- sub("^ecdf_fn_(.*)\\.rds$", "\\1", sa_files)
  sa_names <- intersect(sa_names, valid_sa_names)
  
  individual_banks <- setdiff(sa_names, pooled_banks)
  
  # Optional name filter
  if (!is.null(filter)) {
    individual_banks <- individual_banks[individual_banks %in% filter]
    relevant_pools <- bank_pool_map %>% filter(bank_name %in% filter) %>%
      pull(pool_id) %>% unique()
    pool_representatives <- pool_representatives %>%
      filter(pool_id %in% relevant_pools)
  }
  
  n_individual <- length(individual_banks)
  n_pools      <- nrow(pool_representatives)
  total_markets <- n_individual + n_pools
  
  results_list <- vector("list", total_markets * 2L) # remaining + lost
  j <- 1L
  successful_markets <- 0L
  set.seed(seed)
  
  # ── helper to read ECDF by metric ─────────────────────────────────────────
  get_ecdf <- function(path) {
    obj <- readRDS(path)
    if (is.list(obj)) obj[[metric]] else NULL
  }
  
  # ── loop over INDIVIDUAL (non-pooled) banks ───────────────────────────────
  message("Processing individual (non-pooled) banks …")
  for (bank in individual_banks) {
    
    sa_file   <- file.path(datadir, "SA ECDF Functions",
                           sprintf("ecdf_fn_%s.rds", bank))
    loss_file <- file.path(datadir, "SA ECDF Functions - Loss",
                           sprintf("ecdf_fn_%s.rds", bank))
    
    if (!file.exists(sa_file) || !file.exists(loss_file)) next
    
    sa_ecdf   <- get_ecdf(sa_file)
    loss_ecdf <- get_ecdf(loss_file)
    if (is.null(sa_ecdf) || is.null(loss_ecdf)) next
    
    sa_data   <- knots(sa_ecdf)
    loss_data <- knots(loss_ecdf)
    if (!length(sa_data) || !length(loss_data)) next
    
    sa_n_obs   <- environment(sa_ecdf)$n   %||% length(sa_data)
    loss_n_obs <- environment(loss_ecdf)$n %||% length(loss_data)
    
    ref_mean <- mean(sa_data, na.rm = TRUE)
    
    if (length(sa_data)   > sample_size) sa_data   <- sample(sa_data,   sample_size)
    if (length(loss_data) > sample_size) loss_data <- sample(loss_data, sample_size)
    
    row_weight <- loss_n_obs / length(loss_data)
    
    make_tbl <- function(v, status) tibble(
      market_id       = bank,
      status          = status,
      value           = v,
      normalized_value= v / ref_mean,
      weight          = row_weight
    ) %>% filter(is.finite(normalized_value) & normalized_value > 0)
    
    results_list[[j]] <- make_tbl(sa_data,   "remaining"); j <- j + 1L
    results_list[[j]] <- make_tbl(loss_data, "lost");       j <- j + 1L
    successful_markets <- successful_markets + 1L
  }
  
  # ── loop over POOLED banks ────────────────────────────────────────────────
  message("Processing pooled banks …")
  for (i in seq_len(n_pools)) {
    pool_id           <- pool_representatives$pool_id[i]
    representative_sa <- pool_representatives$bank_name[i]
    
    pooled_sa_file <- file.path(datadir, "Pooled HUC8 ECDFs",
                                sprintf("ecdf_pool_vs_sa_%s.rds", pool_id))
    loss_file      <- file.path(datadir, "SA ECDF Functions - Loss",
                                sprintf("ecdf_fn_%s.rds", representative_sa))
    
    if (!file.exists(pooled_sa_file) || !file.exists(loss_file)) next
    
    pooled_data     <- readRDS(pooled_sa_file)
    sa_ecdf         <- pooled_data$sa_rem_ecdf[[metric]]
    loss_ecdf       <- get_ecdf(loss_file)
    if (is.null(sa_ecdf) || is.null(loss_ecdf)) next
    
    sa_data   <- knots(sa_ecdf)
    loss_data <- knots(loss_ecdf)
    if (!length(sa_data) || !length(loss_data)) next
    
    sa_n_obs   <- environment(sa_ecdf)$n   %||% length(sa_data)
    loss_n_obs <- environment(loss_ecdf)$n %||% length(loss_data)
    
    ref_mean <- mean(sa_data, na.rm = TRUE)
    
    if (length(sa_data)   > sample_size) sa_data   <- sample(sa_data,   sample_size)
    if (length(loss_data) > sample_size) loss_data <- sample(loss_data, sample_size)
    
    row_weight <- loss_n_obs / length(loss_data)
    
    make_tbl <- function(v, status) tibble(
      market_id        = pool_id,
      status           = status,
      value            = v,
      normalized_value = v / ref_mean,
      weight           = row_weight
    ) %>% filter(is.finite(normalized_value) & normalized_value > 0)
    
    results_list[[j]] <- make_tbl(sa_data,   "remaining"); j <- j + 1L
    results_list[[j]] <- make_tbl(loss_data, "lost");       j <- j + 1L
    successful_markets <- successful_markets + 1L
  }
  
  if (successful_markets == 0)
    stop("No markets processed for metric ", metric)
  
  results <- bind_rows(results_list)
  rm(results_list); gc()
  
  # ── full ECDF dataframe (for plotting / stats) ───────────────────────────
  results_ecdf <- results %>%
    group_by(market_id, status) %>%
    arrange(normalized_value) %>%
    mutate(ecdf_y = row_number() / n()) %>%
    ungroup()
  
  loss_mean <- weighted.mean(results$normalized_value[results$status == "lost"],
                             w = results$weight[results$status == "lost"],
                             na.rm = TRUE)
  p_val <- t.test(normalized_value ~ status, data = results)$p.value
  
  list(data = results_ecdf,
       loss_mean = loss_mean,
       p_value = p_val,
       n_markets = successful_markets,
       n_individual = n_individual,
       n_pools = n_pools)
}

# ──────────────────────────────────────────────────────────────────────────────
# 2)  Generic density-PDF plotting function
# ──────────────────────────────────────────────────────────────────────────────
create_wetland_pdf_plot <- function(results_data, loss_mean, p_val,
                                    title        = NULL,
                                    colours      = c(lost = "red",
                                                     remaining = "blue"),
                                    bracket_pad   = 0.10,
                                    star_pad_frac = 0.05,
                                    star_nudge    = -0.15) {  
  sig_stars <- dplyr::case_when(
    p_val < 0.001 ~ "***",
    p_val < 0.01  ~ "**",
    p_val < 0.05  ~ "*",
    TRUE          ~ "ns"
  )
  
  dens_rem  <- density(results_data$normalized_value[results_data$status == "remaining"])
  dens_lost <- density(results_data$normalized_value[results_data$status == "lost"])
  y_max     <- max(dens_rem$y, dens_lost$y, na.rm = TRUE)
  
  y_bracket <- y_max * (1 + bracket_pad)
  y_tick    <- y_max * (1 + bracket_pad - 0.03)
  
  x_left  <- min(1, loss_mean)
  x_right <- max(1, loss_mean)
  
  ggplot(results_data,
         aes(x = normalized_value, colour = status, weight = weight)) +
    stat_density(geom = "line", position = "identity", linewidth = 1.2) +
    geom_vline(xintercept = 1,         linetype = "dashed",
               colour = colours["remaining"], linewidth = .8) +
    geom_vline(xintercept = loss_mean, linetype = "dashed",
               colour = colours["lost"],      linewidth = .8) +
    annotate("segment", x = x_left,  xend = x_right,
             y = y_bracket, yend = y_bracket) +
    annotate("segment", x = x_left,  xend = x_left,
             y = y_bracket, yend = y_tick) +
    annotate("segment", x = x_right, xend = x_right,
             y = y_bracket, yend = y_tick) +
    annotate("text",
             x = (x_left + x_right) / 2 + star_nudge,   
             y = y_bracket + star_pad_frac * y_max,
             label = sig_stars, vjust = 0) +
    scale_x_log10(
      limits = c(0.01, 100),
      breaks = c(0.1, 1, 2, 10, 100, 1000),            
      labels = function(x) format(x, digits = 2, trim = TRUE),
      expand = c(0, 0)
    ) +
    expand_limits(y = y_bracket + (star_pad_frac + 0.01) * y_max) +
    labs(title = title,
         x = NULL,                                     # suppress internal x title
         y = "Density") +
    scale_color_manual(
      values = colours,
      labels = c(lost = "Wetland Loss (1985–2021)",
                 remaining = "Remaining Wetland (2021)"),
      name   = ""
    ) +
    theme_bw(base_size = 14) +
    theme(
      legend.position = "bottom",
      legend.title    = element_text(size = 14),
      legend.text     = element_text(size = 12),
      axis.title      = element_text(size = 16),
      axis.text       = element_text(size = 14),
      plot.title      = element_text(size = 16, hjust = 0.5),
      panel.grid.major = element_blank(),
      panel.grid.minor = element_blank()
    )
}

plot_value <- create_wetland_pdf_plot(res_value$data,
                                      loss_mean = res_value$loss_mean,
                                      p_val     = res_value$p_value,
                                      title     = "Housing Value",
                                      colours   = pal_value,
                                      bracket_pad   = 0.25,
                                      star_pad_frac = 0.12,
                                      star_nudge    = -0.10)

plot_pop  <- create_wetland_pdf_plot(res_pop$data,
                                     loss_mean = res_pop$loss_mean,
                                     p_val     = res_pop$p_value,
                                     title     = "Population",
                                     colours   = pal_pop,
                                     bracket_pad   = 0.25,
                                     star_pad_frac = 0.12,
                                     star_nudge    = -0.10) +
  labs(y = NULL) +                          # remove y title
  theme(axis.title.y = element_blank())

plot_pop
library(patchwork)

pdf_panels <- (plot_value | plot_pop) +
  plot_annotation(
    caption = "Flood Protection Index (Relative to Remaining Service Area Mean)",
    theme = theme(
      plot.caption = element_text(size = 16, hjust = 0.5, vjust = -0.5),
      plot.margin  = margin(b = 20)
    )
  )
pdf_panels

ggsave("CONUS/Graphics/pdfplot_hvpop.png",
       pdf_panels, width = 12, height = 6, dpi = 300)


############################################
############ Land value PDF Plot ###########
############################################


suppressPackageStartupMessages({
  library(dplyr)
  library(purrr)
  library(terra)
  library(sf)
  library(readr)
})

# ------------------------------------------------------------------------------
# Inputs 
# ------------------------------------------------------------------------------

root_dir        <- "C:\Users\indumati\Box\Paper2_final"
sa_mosaic_dir   <- file.path(root_dir, "Service Area Mosaics")
loss_mosaic_dir <- file.path(root_dir, "Service Area Mosaics", "Loss Mosaics")

nolte_landval <- rast(file.path(root_dir, "Land Values", "1 estimates", "places_fmv_vacant.tif"))

sas   <- readRDS(file.path(root_dir, "Service Areas", "ServiceAreas_agg.rds")) |> st_as_sf()
banks <- readRDS(file.path(root_dir, "Footprints", "footprints_and_buffers.rds")) |> st_as_sf()

bank_pool_map <- readRDS(file.path(root_dir, "Extractions and Summaries", "ECDF Functions", "bank_pool_map.rds"))

# ------------------------------------------------------------------------------
# Small utility functions
# ------------------------------------------------------------------------------

norm_name <- function(x) {
  x <- trimws(x)
  x <- tolower(x)
  gsub("_+", "_", x, perl = TRUE)
}

sa_path_for   <- function(id) file.path(sa_mosaic_dir,   paste0(id, "_mosaic.tif"))
loss_path_for <- function(id) file.path(loss_mosaic_dir, paste0(id, "_lossmosaic.tif"))

has_sa_files <- function(id) file.exists(sa_path_for(id)) && file.exists(loss_path_for(id))

summarise_ln <- function(v_ln) {
  v_ln <- v_ln[is.finite(v_ln)]
  tibble(
    n_cells       = length(v_ln),
    mean_ln       = if (length(v_ln)) mean(v_ln) else NA_real_,
    median_ln     = if (length(v_ln)) median(v_ln) else NA_real_,
    mean_usd_ha   = if (length(v_ln)) mean(exp(v_ln)) else NA_real_,
    median_usd_ha = if (length(v_ln)) median(exp(v_ln)) else NA_real_
  )
}

# ------------------------------------------------------------------------------
# Raster alignment (memoized per target raster)
# ------------------------------------------------------------------------------

align_nolte_to <- function(nolte_r, target_r) {
  stopifnot(inherits(nolte_r, "SpatRaster"), inherits(target_r, "SpatRaster"))
  if (!terra::same.crs(nolte_r, target_r)) nolte_r <- terra::project(nolte_r, target_r)
  terra::resample(nolte_r, target_r, method = "bilinear")
}

# Cache aligned Nolte raster per file-backed target raster (by filename + layer)
.nolte_cache <- new.env(parent = emptyenv())

get_nolte_on_grid <- function(nolte_r, mask_r, cache_key) {
  if (exists(cache_key, envir = .nolte_cache, inherits = FALSE)) {
    return(get(cache_key, envir = .nolte_cache, inherits = FALSE))
  }
  aligned <- align_nolte_to(nolte_r, mask_r)
  assign(cache_key, aligned, envir = .nolte_cache)
  aligned
}

# ------------------------------------------------------------------------------
# Extract Nolte land values (ln $/ha) from a mosaic, inside poly, excluding exclude
# ------------------------------------------------------------------------------

extract_landvals_from_mosaic <- function(nolte_r, mosaic_r, poly_sf = NULL,
                                         exclude_sf = NULL, mosaic_layer = 1,
                                         cache_key = NULL) {
  stopifnot(inherits(nolte_r, "SpatRaster"), inherits(mosaic_r, "SpatRaster"))
  
  mask_r <- mosaic_r[[mosaic_layer]]
  
  if (is.null(cache_key)) {
    # fallback cache key if caller doesn't provide one
    cache_key <- paste0("nolte_", terra::ncol(mask_r), "x", terra::nrow(mask_r), "_",
                        paste(terra::res(mask_r), collapse = "x"), "_",
                        terra::crs(mask_r))
  }
  
  nolte_on_grid <- get_nolte_on_grid(nolte_r, mask_r, cache_key)
  
  if (!is.null(poly_sf)) {
    poly_v <- if (inherits(poly_sf, "sf")) terra::vect(poly_sf) else poly_sf
    stopifnot(inherits(poly_v, "SpatVector"))
    
    mask_r        <- terra::crop(mask_r, poly_v)
    nolte_on_grid <- terra::crop(nolte_on_grid, poly_v)
    
    mask_r        <- terra::mask(mask_r, poly_v)
    nolte_on_grid <- terra::mask(nolte_on_grid, poly_v)
  }
  
  if (!is.null(exclude_sf)) {
    ex_v <- if (inherits(exclude_sf, "sf")) terra::vect(exclude_sf) else exclude_sf
    stopifnot(inherits(ex_v, "SpatVector"))
    
    ex_r <- terra::rasterize(ex_v, mask_r, field = 1, background = NA)
    mask_r[!is.na(ex_r)]        <- NA
    nolte_on_grid[!is.na(ex_r)] <- NA
  }
  
  nolte_on_grid <- terra::mask(nolte_on_grid, mask_r)
  vals_ln <- terra::values(nolte_on_grid, na.rm = TRUE)
  vals_ln[is.finite(vals_ln)]
}

# ------------------------------------------------------------------------------
# SOLO market summarizer: loss/bank/remaining for one ID
# ------------------------------------------------------------------------------

summarize_market <- function(id, nolte_r, sas, banks, mosaic_layer = 1) {
  sa_poly   <- sas   |> filter(ID == id)
  bank_poly <- banks |> filter(Name == id)
  
  if (nrow(sa_poly) == 0 || nrow(bank_poly) == 0) return(NULL)
  if (!has_sa_files(id)) return(NULL)
  
  sa_mosaic   <- terra::rast(sa_path_for(id))
  loss_mosaic <- terra::rast(loss_path_for(id))
  
  cache_sa   <- paste0("sa_", id, "_layer", mosaic_layer)
  cache_loss <- paste0("loss_", id, "_layer", mosaic_layer)
  
  lv_bank_ln <- extract_landvals_from_mosaic(
    nolte_r, sa_mosaic, poly_sf = bank_poly, exclude_sf = NULL,
    mosaic_layer = mosaic_layer, cache_key = cache_sa
  )
  
  lv_remaining_ln <- extract_landvals_from_mosaic(
    nolte_r, sa_mosaic, poly_sf = sa_poly, exclude_sf = bank_poly,
    mosaic_layer = mosaic_layer, cache_key = cache_sa
  )
  
  lv_loss_ln <- extract_landvals_from_mosaic(
    nolte_r, loss_mosaic, poly_sf = sa_poly, exclude_sf = bank_poly,
    mosaic_layer = mosaic_layer, cache_key = cache_loss
  )
  
  bind_rows(
    summarise_ln(lv_loss_ln)      |> mutate(group = "loss"),
    summarise_ln(lv_bank_ln)      |> mutate(group = "bank"),
    summarise_ln(lv_remaining_ln) |> mutate(group = "remaining")
  ) |> mutate(
    pool_id   = NA_character_,
    rep_bank  = id,
    market_id = id,
    bank_id   = id,
    .before = 1
  )
}

# ------------------------------------------------------------------------------
# POOL summarizer:
# - pool-level: loss + remaining using SA polygon & mosaics from rep_bank
# - bank-level:  rows per bank using same SA mosaic grid
# ------------------------------------------------------------------------------
summarize_pool <- function(pool_id, rep_bank, bank_names,
                           nolte_r, sas, banks,
                           mosaic_layer = 1) {
  
  sa_poly <- sas |> filter(ID == rep_bank)
  if (nrow(sa_poly) == 0) return(NULL)
  if (!has_sa_files(rep_bank)) return(NULL)
  
  sa_mosaic   <- terra::rast(sa_path_for(rep_bank))
  loss_mosaic <- terra::rast(loss_path_for(rep_bank))
  
  cache_sa   <- paste0("sa_", rep_bank, "_layer", mosaic_layer)
  cache_loss <- paste0("loss_", rep_bank, "_layer", mosaic_layer)
  
  bank_polys <- banks |>
    mutate(Name_norm = norm_name(Name)) |>
    filter(Name_norm %in% norm_name(bank_names))
  
  if (nrow(bank_polys) == 0) return(NULL)
  
  bank_union <- bank_polys |> st_union() |> st_as_sf()
  
  # 1) UNIONED BANKS (on SA mosaic)
  lv_bank_union_ln <- extract_landvals_from_mosaic(
    nolte_r, sa_mosaic, poly_sf = bank_union, exclude_sf = NULL,
    mosaic_layer = mosaic_layer, cache_key = cache_sa
  )
  
  # 2) REMAINING (SA minus unioned banks)
  lv_remaining_ln <- extract_landvals_from_mosaic(
    nolte_r, sa_mosaic, poly_sf = sa_poly, exclude_sf = bank_union,
    mosaic_layer = mosaic_layer, cache_key = cache_sa
  )
  
  # 3) LOSS (loss mosaic within SA, excluding unioned banks)
  lv_loss_ln <- extract_landvals_from_mosaic(
    nolte_r, loss_mosaic, poly_sf = sa_poly, exclude_sf = bank_union,
    mosaic_layer = mosaic_layer, cache_key = cache_loss
  )
  
  bind_rows(
    summarise_ln(lv_bank_union_ln) |> mutate(group = "bank"),
    summarise_ln(lv_loss_ln)       |> mutate(group = "loss"),
    summarise_ln(lv_remaining_ln)  |> mutate(group = "remaining")
  ) |>
    mutate(
      pool_id   = pool_id,
      rep_bank  = rep_bank,
      market_id = rep_bank,
      bank_id   = rep_bank,
      .before = 1
    )
}

# ------------------------------------------------------------------------------
# Safe wrappers
# ------------------------------------------------------------------------------
safe_market <- safely(function(id) summarize_market(id, nolte_landval, sas, banks),
                      otherwise = NULL)
safe_pool <- safely(function(pool_id, rep_bank, bank_names) {
  summarize_pool(pool_id, rep_bank, bank_names, nolte_landval, sas, banks)
}, otherwise = NULL)


# ------------------------------------------------------------------------------
# Build pool membership + solo list 
# ------------------------------------------------------------------------------
pool_lookup <- bank_pool_map |>
  transmute(pool_id = pool_id, bank_name = bank_name) |>
  distinct()

pool_members <- pool_lookup |>
  group_by(pool_id) |>
  summarize(bank_names = list(unique(bank_name)), .groups = "drop") |>
  mutate(
    rep_bank = map_chr(bank_names, ~{
      ok <- .x[map_lgl(.x, has_sa_files)]
      if (length(ok) == 0) NA_character_ else ok[1]
    })
  ) |>
  filter(!is.na(rep_bank))

message("Pools with rep mosaics: ", nrow(pool_members))

all_bank_ids    <- unique(banks$Name)
pooled_bank_ids <- unique(pool_lookup$bank_name)

solo_bank_ids <- setdiff(all_bank_ids, pooled_bank_ids)
solo_bank_ids <- solo_bank_ids[map_lgl(solo_bank_ids, has_sa_files)]
message("Solo banks with mosaics: ", length(solo_bank_ids))

# ------------------------------------------------------------------------------
# Run ALL
# ------------------------------------------------------------------------------
# test = pool_members[1,]
pool_res <- map(seq_len(nrow(pool_members)), function(i) {
  pid <- pool_members$pool_id[i]
  rep <- pool_members$rep_bank[i]
  bns <- pool_members$bank_names[[i]]
  message("[POOL] ", pid, " | rep=", rep, " | n_banks=", length(bns))
  safe_pool(pid, rep, bns)
})
pool_out <- map_dfr(pool_res, "result")

#saveRDS(pool_out,"L:\\Wetland Flood Mitigation\\CONUS\\Land Values\\pooled_landvalue.rds")
solo_res <- map(solo_bank_ids, function(id) {
  message("[SOLO] ", id)
  safe_market(id)
})
solo_out <- map_dfr(solo_res, "result")
#saveRDS(solo_out,"L:\\Wetland Flood Mitigation\\CONUS\\Land Values\\solo_landvalue.rds")
all_out <- bind_rows(pool_out, solo_out)
saveRDS(all_out,"L:\\Wetland Flood Mitigation\\CONUS\\Land Values\\all_landvalue.rds")

message("DONE. Rows: ", nrow(all_out))
print(count(all_out, group))

# ------------------------------------------------------------------------------
# Save outputs + basic failure logs
# ------------------------------------------------------------------------------
out_rds <- file.path(root_dir, "CONUS", "Land Values", "market_landvals_summary_poolaware.rds")
out_csv <- file.path(root_dir, "CONUS", "Land Values", "market_landvals_summary_poolaware.csv")

saveRDS(all_out, out_rds)
write_csv(all_out, out_csv)

pool_fail <- tibble(
  pool_id  = pool_members$pool_id,
  rep_bank = pool_members$rep_bank,
  ok       = map_lgl(pool_res, ~ !is.null(.x$result))
) |> filter(!ok)

solo_fail <- tibble(
  market_id = solo_bank_ids,
  ok        = map_lgl(solo_res, ~ !is.null(.x$result))
) |> filter(!ok)

fail_rds <- file.path(root_dir, "CONUS", "Land Values", "market_landvals_failures_poolaware.rds")
saveRDS(list(pool_fail = pool_fail, solo_fail = solo_fail), fail_rds)

message("Failures | pools: ", nrow(pool_fail), " | solo: ", nrow(solo_fail))


####
# Correcting the bank name mismatches and re-running those
suppressPackageStartupMessages({
  library(dplyr)
  library(purrr)
  library(tibble)
})
ids_rerun <- in_res_not_in_solo

#ids_rerun <- in_res_not_in_solo

# run safe_market again (or reuse rerun_res if you still have it)
rerun_res <- map(ids_rerun, ~ safe_market(.x))

rerun_diag <- tibble(
  id = ids_rerun,
  ok = map_lgl(rerun_res, ~ !is.null(.x$result)),
  err = map_chr(rerun_res, ~ if (!is.null(.x$error)) conditionMessage(.x$error) else NA_character_)) 
%>%
  mutate(
    has_files = map_lgl(id, has_sa_files),
    has_sa    = id %in% sas$ID,
    has_bank  = id %in% banks$Name,
    reason = case_when(
      ok ~ "ok",
      !is.na(err) ~ paste0("ERROR: ", err),
      !has_files ~ "return(NULL): missing mosaic files",
      !has_sa ~ "return(NULL): sas$ID not found",
      !has_bank ~ "return(NULL): banks$Name not found",
      TRUE ~ "return(NULL): unknown (passed simple checks)"
    )
  
)

manual_map <- c(
  "FP_AndL_Everglades_Phase_II_MB" = "FP_AndL_Everglades_Phase_II_MB",
  "MT-Lower_Middle_Yellowstone_Umbrella_Bank-_Ringling_Ranch_Bank_Site" = "MT-Lower/Middle_Yellowstone_Umbrella_Bank-_Ringling_Ranch_Bank_Site",
  "MT-Lower_Middle_Yellowstone_Umbrella-Schwarzkoph-Waage_Bank_Site" = "MT-Lower/Middle_Yellowstone_Umbrella-Schwarzkoph-Waage_Bank_Site",
  "MN_-_Beltrami_-_S06985-1583_Stelter" = "MN_-_Beltrami_-_S06985-1583__Stelter",
  "ND-NDDOT_Kirkeby_Schuster" = "ND-NDDOT_Kirkeby/Schuster",
  "MN_-_Lake_-_S4545-1532_Ziemet_Peterson_Preserve" = "MN_-_Lake_-_S4545-1532_Ziemet/Peterson_Preserve",
  "MN_-_Hennepin_-_S4168-1410_DeLuca_Heritage_Woods" = "MN_-_Hennepin_-_S4168-1410__DeLuca_Heritage_Woods",
  "MN_-_Carver_-_S132-1375_Richards" = "MN_-_Carver_-_S132-1375__Richards",
  "Swallow_Tail_LLC_Nishnabotna_Platte_Rivers_Umbrella_Bank_Site_2" = "Swallow_Tail_LLC_Nishnabotna/Platte_Rivers_Umbrella_Bank_Site_2",
  "Swallow_Tail_LLC_Nishnabotna_Platte_Rivers_Umbrella_Bank_Site_1" = "Swallow_Tail_LLC_Nishnabotna/Platte_Rivers_Umbrella_Bank_Site_1",
  "MN_-_Anoka_-_S3049-1350_Strandlund_Boettcher_Sec_29" = "MN_-_Anoka_-_S3049-1350__Strandlund_Boettcher_Sec_29",
  "UDOT_Northern_Utah_County" = "UDOT_Northern_Utah_County",
  "ND-DU_UMB_Starkweather_Site" = "ND-DU_UMB__Starkweather_Site",
  "Peace_River_Hardee_Co_MB" = "Peace_River/Hardee_Co_MB",
  "Lake_Louisa_Green_Swamp_MB" = "Lake_Louisa/Green_Swamp_MB",
  "Swallow_Tail_LLC_Blackwater_Lamine_Rivers_Umbrella_Bank_Site_1" = "Swallow_Tail_LLC_Blackwater/Lamine_Rivers_Umbrella_Bank_Site_1",
  "Swallow_Tail_LLC_Blackwater_Lamine_Rivers_Umbrella_Bank_Site_2" = "Swallow_Tail_LLC_Blackwater/Lamine_Rivers_Umbrella_Bank_Site_2",
  "ESS_GREEN_1_LLC_Blackwater_Lamine_Rivers_Umbrella_Bank_Site_1" = "ESS_GREEN_1_LLC_Blackwater/Lamine_Rivers_Umbrella_Bank_Site_1",
  "Scuppernong_River_Corridor_Mitigation_Bank" = "Scuppernong_River_Corridor_Mitigation_Bank",
  "Upper_Coastal_Citrus_Co_MB" = "Upper_Coastal/Citrus_Co_MB",
  "WI_-_Monroe_-_H_G_Randall" = "WI_-_Monroe_-_HG_Randall",
  "NE-Fairmont_Rainwater_Basin" = "NE-Fairmont/Rainwater_Basin",
  "Colbert_Cameron_MB" = "Colbert/Cameron_MB",
  "Mitico_LLC_Grand_Chariton_EDU_Umbrella_Bank_Site_1" = "Mitico_LLC_Grand/Chariton_EDU_Umbrella_Bank_Site_1",
  "PSUMBI_-_East_Branch_Codorus_Creek_Bank" = "PSUMBI_-_East_Branch_Codorus_Creek_Bank",
  "Swallow_Tail_LLC_Grand_Chariton_Rivers_Umbrella_Bank_Site_1" = "Swallow_Tail_LLC_Grand/Chariton_Rivers_Umbrella_Bank_Site_1",
  "R_A_Burgess_Mitigation_Bank" = "RA_Burgess_Mitigation_Bank",
  "Barra_Farms-Phase_I" = "Barra_Farms-Phase_I",
  "Swallow_Tail_LLC_Osage_South_Grand_Rivers_Wetland_and_Stream_Umbrella_Mitigation_Bank_Site_One" = "Swallow_Tail_LLC_Osage/South_Grand_Rivers_Wetland_and_Stream_Umbrella_Mitigation_Bank_Site_One",
  "Lick_Creek_Wetland_Bank_2" = "Lick_Creek_Wetland_Bank__2",
  "Lick_Creek_Wetland_Bank_1" = "Lick_Creek_Wetland_Bank__1",
  "West_Kentucky_Wetwoods_Mitigation_Bank_2" = "West_Kentucky_Wetwoods_Mitigation_Bank_2",
  "Beverly_Mitigation_Bank_Site_1" = "Beverly_Mitigation_Bank_Site__1",
  "Honey_Springs_Mitigation_Bank_-_Elk_Tract_Old_Field_Tract" = "Honey_Springs_Mitigation_Bank_-_Elk_Tract/Old_Field_Tract",
  "Buttahatchie_River_Mississippi_Phase_I" = "Buttahatchie_River/Mississippi_Phase_I"
)

banks <- banks %>% mutate(Name_orig = Name)

inv_map <- setNames(names(manual_map), unname(manual_map))  # current -> desired
banks$Name[banks$Name %in% names(inv_map)] <- inv_map[banks$Name[banks$Name %in% names(inv_map)]]

# ------------------------------------------------------------------------------
# 2) Rerun solo_out for markets that are in res but missing from solo_out (exclude pools)
# ------------------------------------------------------------------------------
ids_to_rerun <- in_res_not_in_solo[!grepl("^Pool_", in_res_not_in_solo)]
ids_to_rerun <- ids_to_rerun[map_lgl(ids_to_rerun, has_sa_files)]

message("Rerunning ", length(ids_to_rerun), " solo markets...")

rerun_res <- map(ids_to_rerun, function(id) {
  message("[RERUN SOLO] ", id)
  safe_market(id)
})

solo_rerun_out <- map_dfr(rerun_res, "result")



# ------------------------------------------------------------------------------
# 3) Append + remove duplicates (keep one row per bank_id × group)
# ------------------------------------------------------------------------------
solo_out <- bind_rows(solo_out, solo_rerun_out) %>%
  arrange(bank_id, group) %>%
  distinct(bank_id, group, .keep_all = TRUE)

# sanity check: should be no duplicates left
solo_out %>% count(bank_id, group) %>% filter(n > 1)

### Make plot ---------------------------------------------------

create_land_value_plot_with_tests <- function(all_out,
                                              bracket_pad = 0.30,
                                              star_pad_frac = 0.10) {
  suppressPackageStartupMessages({
    library(dplyr)
  })
  
  # Prepare normalized data (only remaining and loss)
  remaining_means <- all_out %>%
    filter(group == "remaining") %>%
    select(rep_bank, mean_usd_ha) %>%
    rename(remaining_mean = mean_usd_ha)
  
  plot_data <- all_out %>%
    filter(group %in% c("remaining", "loss")) %>%
    left_join(remaining_means, by = "rep_bank") %>%
    filter(!is.na(mean_usd_ha), !is.na(remaining_mean), remaining_mean > 0) %>%
    mutate(
      normalized_value = mean_usd_ha / remaining_mean,
      group_label = case_when(
        group == "remaining" ~ "Remaining Wetland (2021)",
        group == "loss"      ~ "Wetland Loss (1985–2021)"
      )
    ) %>%
    filter(is.finite(normalized_value) & normalized_value > 0)
  
  # Summary stats
  group_stats <- plot_data %>%
    group_by(group_label) %>%
    summarise(
      mean_val   = mean(normalized_value, na.rm = TRUE),
      median_val = median(normalized_value, na.rm = TRUE),
      n          = n(),
      .groups    = "drop"
    )
  
  # T-test (unweighted)
  p_val <- t.test(normalized_value ~ group, data = plot_data)$p.value
  
  sig_stars <- case_when(
    p_val < 0.001 ~ "***",
    p_val < 0.01  ~ "**",
    p_val < 0.05  ~ "*",
    TRUE          ~ "ns"
  )
  
  # Create plot using your fixed plot function
  p <- create_land_value_distribution_plot(
    all_out,
    bracket_pad   = bracket_pad,
    star_pad_frac = star_pad_frac
  )
  
  list(
    plot         = p,
    statistics   = group_stats,
    p_value      = p_val,
    significance = sig_stars,
    data         = plot_data
  )
}

# Usage:
result <- create_land_value_plot_with_tests(all_out)
print(result$plot)
# result$statistics

ggsave("L:/Wetland Flood Mitigation/CONUS/Graphics/pdfplot_landvalue.png", result$plot, width=12, height=8, dpi=300)
result$statistics



#############################
######## FIGURE TOP 10 ######
#############################

library(tidyverse)
sa_dir <- "Extractions and Summaries/Extract Service Areas"
bank_dir <- "Extractions and Summaries/Extract Banks"
loss_dir <- "Extractions and Summaries/Extract Loss"
# ----------------------------
# Define top 10 bank name groups
# ----------------------------
top10_names <- list(
  c("Loxahatchee_MB"),
  c("Little_Pine_Island_MB"),
  c("Panther_Island_MB", "Panther_Island_MB_-_Expansion"),
  c("Florida_Wetlandsbank_at_Pembroke_Pines_MB"),
  c("Big_Cypress_MB_Phase_I-V", "Big_Cypress_MB_Phase_VI"),
  c("Crosby_Island_Marsh_MB"),
  c("St_Johns_MB"),
  c("Katy_Prairie_Stream"),
  c("Poa_Boy_MB"),
  c("FP_AndL_Everglades_Phase_I_MB", "FP_AndL_Everglades_Phase_II_MB")
)

# ----------------------------
# Updated load helper (vectorized match)
# ----------------------------
load_and_rename <- function(dir, prefix, names_vec) {
  files <- list.files(dir, pattern = paste0("^", prefix, ".*\\.rds$"), full.names = TRUE)
  name_pattern <- paste0(names_vec, collapse = "|")
  
  files |> 
    keep(~ str_detect(.x, name_pattern)) |> 
    map(readRDS) |> 
    map(rename_cols)
}

# ----------------------------
# Updated file path matcher
# ----------------------------
get_matching_files <- function(dir, prefix, names_vec) {
  files <- list.files(dir, pattern = paste0("^", prefix, ".*\\.rds$"), full.names = TRUE)
  name_pattern <- paste0(names_vec, collapse = "|")
  keep(files, ~ str_detect(.x, name_pattern))
}

# ----------------------------
# Load and rename all RDS files
# ----------------------------
sa_rds_list   <- map(top10_names, ~ load_and_rename(sa_dir,   "sa_extract_",   .x))
bank_rds_list <- map(top10_names, ~ load_and_rename(bank_dir, "bank_extract_", .x))
loss_rds_list <- map(top10_names, ~ load_and_rename(loss_dir, "sa_extract_",   .x))

# ----------------------------
# Get all relevant file paths
# ----------------------------
sa_files   <- map(top10_names, ~ get_matching_files(sa_dir,   "sa_extract_",   .x)) |> flatten_chr()
bank_files <- map(top10_names, ~ get_matching_files(bank_dir, "bank_extract_", .x)) |> flatten_chr()
loss_files <- map(top10_names, ~ get_matching_files(loss_dir, "sa_extract_",   .x)) |> flatten_chr()

# ----------------------------
# Build a file metadata table
# ----------------------------
files_df <- tibble(
  type = c(rep("bank", length(bank_files)),
           rep("service_area", length(sa_files)),
           rep("loss", length(loss_files))),
  file = c(bank_files, sa_files, loss_files)
) |> 
  mutate(
    full_name = basename(file) |> str_remove("\\.rds$"),
    base = full_name |> 
      str_remove("^(bank|sa)_extract_") |>
      str_replace("_-.*", "") |>
      str_replace("_Phase.*", "") |>
      str_replace("_I-VI?$", "")
  )
# ----------------------------
# Combine multiple bank files per base
# ----------------------------
bank_data <- files_df |> 
  filter(type == "bank") |> 
  group_by(base) |> 
  summarise(
    canonical = first(full_name),
    bank = list(map(file, ~ readRDS(.x) |> rename_cols()) |> bind_rows()),
    .groups = "drop"
  )

# ----------------------------
# Use first matching service area per base
# ----------------------------
sa_data <- files_df |>
  filter(type == "service_area") |>
  group_by(base) |>
  slice(1) |>
  mutate(service_area = map(file, ~ readRDS(.x) |> rename_cols())) |>
  dplyr::select(base, service_area)

# ----------------------------
# Use first matching loss per base
# ----------------------------
loss_data <- files_df |>
  filter(type == "loss") |>
  group_by(base) |>
  slice(1) |>
  mutate(loss = map(file, ~ readRDS(.x) |> rename_cols())) |>
  dplyr::select(base, loss)

# ----------------------------
# Merge all data into one nested table
# ----------------------------
data_nested <- bank_data |>
  left_join(sa_data, by = "base") |>
  left_join(loss_data, by = "base")

# ----------------------------
# Build long-form tidy data for plotting
# ----------------------------
plot_df <- data_nested |> 
  mutate(plot = pmap(
    list(canonical, bank, service_area, loss),
    \(canon, bank_df, sa_df, loss_df) {
      bind_rows(
        loss_df |> mutate(group = "loss"),
        sa_df   |> mutate(group = "service_area"),
        bank_df |> mutate(group = "bank")
      ) |>
        mutate(
          bank_name = canon |> 
            str_remove("^bank_extract_") |> 
            str_replace_all("_", " "),
          log_hu = log10(housing_units + 1)
        )
    }
  )) |> 
  pull(plot) |> 
  bind_rows()

# ----------------------------
# sample large groups for plotting
# ----------------------------
max_rows <- 10000

plot_df_sampled <- plot_df %>%
  group_by(group) %>%
  mutate(n_rows = n()) %>%
  group_modify(~ if (.x$n_rows[1] > max_rows) slice_sample(.x, n = max_rows) else .x) %>%
  ungroup() %>%
  dplyr::select(-n_rows)


bank_name_labels <- c(
  "Loxahatchee MB" = "Loxahatchee",
  "Little Pine Island MB" = "Little Pine Island",
  "Panther Island MB" = "Panther Island",
  "Florida Wetlandsbank at Pembroke Pines MB" = "Pembroke Pines",
  "Big Cypress MB Phase I-V" = "Big Cypress",
  "Crosby Island Marsh MB" = "Crosby Island",
  "St Johns MB" = "St. Johns",
  "Katy Prairie Stream" = "Katy Prairie",
  "Poa Boy MB" = "Poa Boy",
  "FP AndL Everglades Phase I MB" = "Everglades")

library(ggplot2)

ggplot(plot_df, aes(x = log_hu, colour = group, fill = group)) +
  geom_density(alpha = 0.25) +
  facet_wrap(~ bank_name, scales = "free", ncol = 3) +
  labs(
    x = "Log downstream housing units",
    y = "Density"
  ) +
  theme_classic()

# Nicer plot
ggplot(plot_df, aes(x = log_hu, colour = group, fill = group)) +
  stat_density(geom = "line", position = "identity", linewidth = 1.2) +
  facet_wrap(~ bank_name, scales = "free", ncol = 3,
             labeller = labeller(bank_name = bank_name_labels)) +
  labs(
    x = "Log downstream housing units",
    y = "Density"
  ) +
  scale_color_manual(
    values = c(
      loss = "red",
      service_area = "blue",
      bank = "#00BA38"
    ),
    labels = c(
      loss = "Wetland Loss (1985–2021)",
      service_area = "Remaining Wetland (2021)",
      bank = "Bank (2021)"
    ),
    name = ""
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


##########################
######## FIGURE S3 #######
##########################

ecdf_loss = readRDS(paste0(datadir, "Extractions and Summaries/ECDF Functions/Loss ECDF Functions/ecdf_fn_Big_Cypress_MB_Phase_I-V.rds"))
ecdf_sa = readRDS(paste0(datadir, "Extractions and Summaries/ECDF Functions/SA ECDF Functions/ecdf_fn_Big_Cypress_MB_Phase_I-V.rds"))
ecdf_bank = readRDS("Extractions and Summaries/ECDF Functions/Bank ECDF Functions/bank_ecdf_fn_Big_Cypress_MB_Phase_I-V.rds")


# Helper: Safely convert ECDF to data frame if it's valid
ecdf_to_df_safe <- function(ecdf_obj, varname, source_label) {
  if (is.function(ecdf_obj)) {
    x_vals <- environment(ecdf_obj)$x
    y_vals <- sapply(x_vals, ecdf_obj)
    return(data.frame(
      value = x_vals,
      ecdf = y_vals,
      variable = varname,
      source = source_label
    ))
  } else {
    return(NULL)
  }
}

# Combine only valid ECDFs
df_all <- bind_rows(
  ecdf_to_df_safe(ecdf_bank$housing_units, "Housing Units", "Bank"),
  ecdf_to_df_safe(ecdf_loss$housing_units, "Housing Units", "Loss"),
  ecdf_to_df_safe(ecdf_sa$housing_units, "Housing Units", "Service Area"),
  
  ecdf_to_df_safe(ecdf_bank$housing_value, "Housing Value", "Bank"),
  ecdf_to_df_safe(ecdf_loss$housing_value, "Housing Value", "Loss"),
  ecdf_to_df_safe(ecdf_sa$housing_value, "Housing Value", "Service Area"),
  
  ecdf_to_df_safe(ecdf_bank$population, "Population", "Bank"),
  ecdf_to_df_safe(ecdf_loss$population, "Population", "Loss"),
  ecdf_to_df_safe(ecdf_sa$population, "Population", "Service Area")
)

# Clean labels
df_all$variable <- factor(df_all$variable,
                          levels = c("Housing Units", "Housing Value", "Population"))


# HV ---

df_hv = df_all %>% 
  filter(variable=="Housing Value")

ecdfplot_hv = ggplot(df_hv, aes(x = value, y = ecdf, color = source)) +
  geom_line(linewidth = 1.2) +
  scale_x_continuous(
    labels = label_dollar(scale = 1e-9, suffix = "B"),
    expand = c(0, 0)
  )+
  scale_color_manual(
    values = c("Loss" = "red", "Service Area" = "blue", "Bank" = "#00BA38"),
    labels = c(
      "Bank" = "Bank",
      "Service Area" = "Remaining Service Area (2021)",
      "Loss" = "Service Area Loss (1985–2021)"
    ),
    name = ""
  ) +
  labs(
    title = "Big Cypress Wetland Mitigation Bank",
    x = "Downstream Total Housing Value",
    y = "Cumulative Probability"
  ) +
  theme_bw(base_size = 14) +
  theme(
    legend.position = "bottom",
    legend.text = element_text(size = 20),
    axis.title = element_text(size = 20),
    axis.text = element_text(size = 20),
    plot.title = element_text(size = 25, hjust = 0.5),
    panel.grid.major = element_blank(),
    panel.grid.minor = element_blank(),
    plot.margin = margin(t = 10, r = 50, b = 10, l = 10))  # increase left/right margins

ecdfplot_hv


# POP ---

df_pop = df_all %>% 
  filter(variable=="Population")

ecdfplot_pop = ggplot(df_pop, aes(x = value, y = ecdf, color = source)) +
  geom_line(linewidth = 1.2) +
  scale_x_continuous(
    labels = comma,
    expand = c(0,0)
  ) +
  scale_color_manual(
    values = c("Loss" = "red", "Service Area" = "blue", "Bank" = "#00BA38"),
    labels = c(
      "Bank" = "Bank",
      "Service Area" = "Remaining Service Area (2021)",
      "Loss" = "Service Area Loss (1985–2021)"
    ),
    name = ""
  ) +
  labs(
    title = "Big Cypress Wetland Mitigation Bank",
    x = "Downstream Population",
    y = "Cumulative Probability"
  ) +
  theme_bw(base_size = 14) +
  theme(
    legend.position = "bottom",
    legend.text = element_text(size = 20),
    axis.title = element_text(size = 20),
    axis.text = element_text(size = 20),
    plot.title = element_text(size = 25, hjust = 0.5),
    panel.grid.major = element_blank(),
    panel.grid.minor = element_blank(),
    plot.margin = margin(t = 10, r = 50, b = 10, l = 10))  # increase left/right margins

ecdfplot_pop

##########################
######## TABLE 1  ########
##########################

summary_loss = readRDS(paste0("Extractions and Summaries/loss_summary.rds"))
summary_bank = readRDS(paste0("Extractions and Summaries/bank_summary.rds"))
sa_summary = readRDS(paste0("Extractions and Summaries/sa_summary.rds"))

summary_bank <- summary_bank %>%
  mutate(bank_name= recode(bank_name, 
                           "FP_AndL_Everglades_Phase_II_MB" = "FP_L_Everglades_Phase_II_MB"))
top10_names <- list(
  "Loxahatchee_MB",
  "Little_Pine_Island_MB",
  "Panther_Island_MB", 
  "Panther_Island_MB_-_Expansion",
  "Florida_Wetlandsbank_at_Pembroke_Pines_MB",
  "Big_Cypress_MB_Phase_I-V", "Big_Cypress_MB_Phase_VI",
  "Crosby_Island_Marsh_MB",
  "St_Johns_MB",
  "Katy_Prairie_Stream",
  "Poa_Boy_MB",
  "FP_L_Everglades_Phase_II_MB"
)

summary_loss_top10  <- summary_loss  %>% filter(sa_id %in% top10_names)
summary_bank_top10  <- summary_bank  %>% filter(bank_name %in% top10_names)
summary_bank_top10 <-summary_bank_top10 %>% 
  rename(sa_id = bank_name)
summary_sa_top10  <- sa_summary  %>% filter(sa_id %in% top10_names)

missing_banks <- setdiff(top10_names, summary_bank_top10$sa_id)
missing_banks

# Add a source column to each summary
summary_loss_top10 <- summary_loss_top10 %>%
  mutate(source = "Loss") %>%
  dplyr::select(sa_id, sum_housing_units, mean_housing_units, max_housing_units, source)

summary_bank_top10 <- summary_bank_top10 %>%
  mutate(source = "Bank") %>%
  dplyr::select(sa_id, sum_housing_units, mean_housing_units, max_housing_units, source)

summary_sa_top10 <- summary_sa_top10 %>%
  mutate(source = "Service Area") %>%
  dplyr::select(sa_id, sum_housing_units, mean_housing_units, max_housing_units, source)

# Combine all
summary_combined <- bind_rows(
  summary_loss_top10,
  summary_bank_top10,
  summary_sa_top10
)

# Sort by descending sum_housing_units
summary_sorted <- summary_combined %>%
  arrange(desc(sum_housing_units))

# Round large numbers for better readability (optional)
summary_sorted_fmt <- summary_sorted %>%
  mutate(across(contains("housing_units"), ~ scales::comma(.)))

# Output to LaTeX
kable(summary_sorted_fmt, format = "latex", booktabs = TRUE,
      caption = "Housing Units Summary for South Florida Wetland Areas",
      col.names = c("Bank / Area", "Sum", "Mean", "Max", "Type")) %>%
  kable_styling(latex_options = c("striped", "hold_position"))

#############################
#### T Test table S2 and S3############
########################################

# Check for autocorr:
# ------------------------------------------------------------------
# libraries  (add SpatialPack and sf)
# ------------------------------------------------------------------
library(dplyr)
library(purrr)
library(broom)
library(glue)
library(sf)          # to drop geometry
library(SpatialPack) # for modified.ttest

# ------------------------------------------------------------------
# helper: read, drop geometry, make column names unique & informative
# ------------------------------------------------------------------
safe_read_rds_renamed <- function(path) {
  df <- tryCatch(readRDS(path), error = function(e) NULL)
  if (is.null(df)) return(NULL)
  
  # remove sf geometry if present
  if (inherits(df, "sf")) df <- st_drop_geometry(df)
  
  # make sure we have at least three metrics to rename
  if (ncol(df) >= 4) {
    names(df)[2:4] <- c("housing_units", "housing_value", "population")
  }
  df
}

# ------------------------------------------------------------------
# main wrapper around SpatialPack::modified.ttest
# ------------------------------------------------------------------
run_modified_t_tests_for_bank <- function(bank_names,
                                          max_pts = 10000,  # target size after sampling
                                          seed    = 1L) {
  
  set.seed(seed)
  canon <- bank_names[1]
  
  # ----- build paths ------------------------------------------------
  bank_paths <- map(bank_names,
                    ~ file.path("Extractions and Summaries", "Extractions with Geom",
                                "Bank Extract with Geom",
                                glue("bank_extract_{.x}withgeom.rds")))
  
  sa_path   <- file.path("Extractions and Summaries", "Extractions with Geom",
                         "SA Extract with Geom",
                         glue("sa_extract_{canon}.rds"))
  
  # loss lives in the same “sa_extract” folder according to your note
  loss_path <- file.path("Extractions and Summaries", "Extractions with Geom",
                         "Loss Extract with Geom",
                         glue("sa_extract_{canon}.rds"))
  
  # ----- read data --------------------------------------------------
  bank_df <- bank_paths %>% map(safe_read_rds_renamed) %>% compact() %>% bind_rows()
  sa_df   <- safe_read_rds_renamed(sa_path)
  loss_df <- safe_read_rds_renamed(loss_path)
  
  # bail out early if something is missing
  if (any(vapply(list(bank_df, sa_df, loss_df), is.null, logical(1)))) return(NULL)
  
  metrics <- c("housing_units", "housing_value", "population")
  
  # helper for one comparison (bank or service-area)
  compare_group <- function(group_df, label) {
    
    map_dfr(metrics, function(metric) {
      
      # skip completely NA metrics
      if (all(is.na(group_df[[metric]]))) return(NULL)
      
      # ----------------------------------------------------------------
      # Assemble a single data.frame on which to run modified.ttest:
      #  • val  – numeric metric
      #  • grp  – binary indicator: 1 = loss, 0 = comparison group
      #  • x,y  – coordinates (already numeric in your example)
      # ----------------------------------------------------------------
      comb <- bind_rows(
        loss_df  %>% select(x, y, val = !!sym(metric)) %>% mutate(grp = 1),
        group_df %>% select(x, y, val = !!sym(metric)) %>% mutate(grp = 0)
      ) |>
        filter(!is.na(val))
      
      # down-sample if the test set is large
      if (nrow(comb) > max_pts) comb <- sample_n(comb, max_pts)
      
      # SpatialPack needs a matrix of coordinates
      coords <- as.matrix(comb[, c("x", "y")])
      
      # ----------------------------------------------------------------
      # modified.ttest: correlation between metric (val) and group dummy
      # The point-biserial r is numerically equivalent to a 2-sample t,
      # so this is the spatially corrected analogue of your original test. :contentReference[oaicite:0]{index=0}
      # ----------------------------------------------------------------
      mod <- SpatialPack::modified.ttest(comb$val, comb$grp, coords, nclass=NULL)
      
      tibble(
        bank       = canon,
        metric     = metric,
        comparison = label,
        mean_loss  = mean(loss_df[[metric]],  na.rm = TRUE),
        mean_group = mean(group_df[[metric]], na.rm = TRUE),
        corr       = mod$corr,
        Fstat      = mod$Fstat,
        dof        = mod$dof,
        p_value    = formatC(mod$p.value, format = "e", digits = 10),
        sampled_n  = nrow(comb)
      )
    })
  }
  
  bind_rows(
    compare_group(bank_df, "bank"),
    compare_group(sa_df,   "service_area")
  )
}

# ------------------------------------------------------------------
# run everything exactly as before
# ------------------------------------------------------------------
top10_names <- list(
  c("Loxahatchee_MB"),
  c("Little_Pine_Island_MB"),
  c("Panther_Island_MB", "Panther_Island_MB_-_Expansion"),
  c("Florida_Wetlandsbank_at_Pembroke_Pines_MB"),
  c("Big_Cypress_MB_Phase_I-V", "Big_Cypress_MB_Phase_VI"),
  c("Crosby_Island_Marsh_MB"),
  c("St_Johns_MB"),
  c("Katy_Prairie_Stream"),
  c("Poa_Boy_MB"),
  c("FP_L_Everglades_Phase_I_MB", "FPL_Everglades_Phase_II_MB"))

# all_modtests <- map_dfr(top10_names, run_modified_t_tests_for_bank)

result <- map_dfr(FP, run_modified_t_tests_for_bank)
result
# Example: housing-unit rows only
result_hu <- result %>% filter(metric == "housing_units")
head(result_hu)


result_hvpop = result %>% filter(metric %in% c("housing_value", "population"))


############# ONE TAILED

safe_read_rds_renamed <- function(path) {
  df <- tryCatch(readRDS(path), error = function(e) NULL)
  if (!is.null(df) && ncol(df) >= 4) {
    colnames(df)[2:4] <- c("housing_units", "housing_value", "population")
  }
  df
}

# ------------------------------------------------------------------
# main wrapper around SpatialPack::modified.ttest with one-tailed p
# ------------------------------------------------------------------
run_modified_t_tests_for_bank <- function(bank_names,
                                          max_pts = 10000,
                                          seed = 1L) {
  
  set.seed(seed)
  canon <- bank_names[1]
  
  bank_paths <- map(bank_names,
                    ~ file.path("Extractions and Summaries", "Extractions with Geom",
                                "Bank Extract with Geom",
                                glue("bank_extract_{.x}withgeom.rds")))
  
  sa_path   <- file.path("Extractions and Summaries", "Extractions with Geom",
                         "SA Extract with Geom",
                         glue("sa_extract_{canon}.rds"))
  
  loss_path <- file.path("Extractions and Summaries", "Extractions with Geom",
                         "Loss Extract with Geom",
                         glue("sa_extract_{canon}.rds"))
  
  bank_df <- bank_paths %>% map(safe_read_rds_renamed) %>% compact() %>% bind_rows()
  sa_df   <- safe_read_rds_renamed(sa_path)
  loss_df <- safe_read_rds_renamed(loss_path)
  
  if (any(vapply(list(bank_df, sa_df, loss_df), is.null, logical(1)))) return(NULL)
  
  metrics <- c("housing_units", "housing_value", "population")
  
  compare_group <- function(group_df, label) {
    map_dfr(metrics, function(metric) {
      if (all(is.na(group_df[[metric]]))) return(NULL)
      
      comb <- bind_rows(
        loss_df  %>% dplyr::select(x, y, val = !!sym(metric)) %>% mutate(grp = 1),
        group_df %>% dplyr::select(x, y, val = !!sym(metric)) %>% mutate(grp = 0)
      ) |>
        filter(!is.na(val))
      
      if (nrow(comb) > max_pts) comb <- sample_n(comb, max_pts)
      coords <- as.matrix(comb[, c("x", "y")])
      
      mod <- SpatialPack::modified.ttest(comb$val, comb$grp, coords, nclass = NULL)
      
      # One-tailed adjustment: loss > comparison (i.e., positive correlation)
      one_tailed_p <- if (mod$corr > 0) mod$p.value / 2 else 1 - (mod$p.value / 2)
      
      tibble(
        bank         = canon,
        metric       = metric,
        comparison   = label,
        mean_loss    = mean(loss_df[[metric]],  na.rm = TRUE),
        mean_group   = mean(group_df[[metric]], na.rm = TRUE),
        corr         = mod$corr,
        Fstat        = mod$Fstat,
        dof          = mod$dof,
        p_value_1tail = formatC(one_tailed_p, format = "e", digits = 10),
        sampled_n    = nrow(comb)
      )
    })
  }
  
  bind_rows(
    compare_group(bank_df, "bank"),
    compare_group(sa_df,   "service_area")
  )
}

# ------------------------------------------------------------------
# Run tests (example with top 10)
# ------------------------------------------------------------------
top10_names <- list(
  c("Loxahatchee_MB"),
  c("Little_Pine_Island_MB"),
  c("Panther_Island_MB", "Panther_Island_MB_-_Expansion"),
  c("Florida_Wetlandsbank_at_Pembroke_Pines_MB"),
  c("Big_Cypress_MB_Phase_I-V", "Big_Cypress_MB_Phase_VI"),
  c("Crosby_Island_Marsh_MB"),
  c("St_Johns_MB"),
  c("Katy_Prairie_Stream"),
  c("Poa_Boy_MB"),
  c("FP_L_Everglades_Phase_I_MB", "FP_L_Everglades_Phase_II_MB"))

# Example: single or multiple bank group
result_onetail <- map_dfr(top10_names, run_modified_t_tests_for_bank)
result_onetail

# Filter examples
result_onetail_hu     <- result_onetail %>% filter(metric == "housing_units")
result_hvpop  <- result %>% filter(metric %in% c("housing_value", "population"))

result_onetail_hv = result_onetail %>% filter(metric == "housing_value")
result_onetail_pop = result_onetail %>% filter(metric == "population")

result_onetail_hv
View(result_onetail_hv)

View(result_onetail_pop)

