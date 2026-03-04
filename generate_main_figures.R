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
})

pth <- function(...) file.path("C:/Users/indumati/Box/Paper2_final", ...)

######################################
######          FIGURE 1       ######
######                          ######
######################################



sanitize_name <- function(name) {
  name <- gsub("[^A-Za-z0-9_-]", "_", name)  # Replace non-alphanumeric characters with "_"
  name <- gsub("_+", "_", name)  # Replace multiple "_" with a single "_"
  name <- gsub("^_|_$", "", name)  # Remove leading or trailing "_"
  return(name)
}

# Load relevant data

usa = st_read(pth("Geo Boundaries", "USA States", "tl_2020_us_state.shp"))
  non_conus <- c("15", "02", "60", "66", "69", "72", "78")
  conus <- usa[!usa$STATEFP %in% non_conus, ]
  conus = st_transform(conus, st_crs(serviceareas))
  conus = vect(conus)
mb = readRDS(pth("Footprints", "mb_pointdata.rds"))
  mb = st_transform(mb, st_crs(serviceareas))
  mb$Name <- sapply(mb$Name, sanitize_name)
  mb = mb %>% 
    filter(Name %in% summary_loss$sa_id)
serviceareas = readRDS(pth("Service Areas", "ServiceAreas_agg.rds"))
  serviceareas = terra::intersect(serviceareas, conus)
  serviceareas = st_as_sf(serviceareas)
  serviceareas = serviceareas %>% 
    filter(ID %in% mb$Name)
summary_loss = readRDS(pth("Extractions and Summaries", "Summaries Loss", "loss_summary.rds"))



# Rename columns for clarity (optional)
colnames(mb)[colnames(mb) == "X"] <- "longitude"
colnames(mb)[colnames(mb) == "Y"] <- "latitude"

# Convert to an sf object with coordinates
mb_sf <- st_as_sf(mb, coords = c("longitude", "latitude"), crs = 4326)
mb_sf = st_transform(mb_sf, st_crs(serviceareas))

mb_sf = st_crop(mb_sf, conus)

mb_sf <- st_centroid(mb)


###

ggplot() +
  # CONUS boundaries
  geom_sf(data = conus, fill = "white", color = "black", size = 0.5, aes(fill = "CONUS")) +
  # Service areas
  geom_sf(data = serviceareas, fill = "gray", color = "gray18", alpha = 0.5, aes(fill = "Service Areas")) + 
  # Mitigation bank points
  geom_sf(data = mb_sf, color = "dodgerblue3", size = 1.5, aes(color = "Mitigation Banks")) +
  
  # North arrow
  annotation_north_arrow(
    location = "bl", which_north = "true",
    style = north_arrow_fancy_orienteering(),
    height = unit(1.5, "cm"), width = unit(1, "cm")
  ) +
  
  # Scale bar (max 1000 km, 4 breaks at 0, 250, 500, 750, 1000)
  annotation_scale(
    location = "bl", width_hint = 0.25,
    bar_cols = c("black", "white"),
    text_cex = 0.8, 
    breaks = c(0, 250, 500, 750, 1000)
  ) +
  
  # Clean theme
  theme_void()

# Inset added in ArcGIS 

##############################################################################
######                                                                ########
######          FIGURE 2 (old version, no pooled service areas)       ########
######                                                                ########
##############################################################################

suppressPackageStartupMessages({
  library(tidyverse)
  library(ggplot2)
  library(scales)
  library(patchwork)
  library(vctrs)
})

compare_remaining_vs_lost <- function(datadir, sample_size = 10000, seed = 123, filter = NULL) {
  suppressPackageStartupMessages({
    library(tidyverse)
    library(scales)
    library(vctrs)
  })
  
  message("Scanning directory for bank ECDF files …")
  sa_files <- list.files(file.path(datadir, "SA ECDF Functions"), full.names = FALSE)
  sa_names <- sub("^ecdf_fn_(.*)\\.rds$", "\\1", sa_files)
  
  
  # ✅ Filter sa_names if filter is provided
  if (!is.null(filter)) {
    sa_names <- sa_names[sa_names %in% filter]
  }
  
  total_banks <- length(sa_names)
  if (total_banks == 0) stop("No ECDF files found for specified banks in directory: ", datadir)
  
  results_list <- vector("list", total_banks * 2)
  j <- 1
  successful_banks <- 0
  set.seed(seed)
  
  for (i in seq_along(sa_names)) {
    bank <- sa_names[i]
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
    
    row_weight <- loss_n_obs / length(loss_data)   # each row represents this many m² lost
    
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
  y_bracket <- y_max * 1.05      # horizontal bar
  y_tick    <- y_max * 0.97      # vertical “ticks” downward
  
  # keep bracket direction left-to-right
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
    ) +     #  ← new
    expand_limits(y = y_bracket + 0.06 * y_max) +
    labs(x = "Value (Relative to Remaining Service Area Mean)",
         y = "Density"
         # title = paste0("Pooled PDFs of Remaining and Lost Wetlands in Mitigation Bank Service Areas\n",
         # "(weighted by Lost Wetland Area)")
    ) +
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
  
  # ecdf_plot <- ggplot(results,
  #                     aes(x = normalized_value, colour = status, weight = weight)) +
  #   stat_ecdf(geom = "step", linewidth = 1.2) +
  #   scale_x_log10(
  #     limits = c(0.01, 100),
  #     breaks = c(0.1, 1, 10, 100),
  #     labels = scales::label_number(accuracy = 1, trim = TRUE),
  #     expand = c(0, 0)
  #   ) +      #  ← same tweak here
  #   labs(x = "Value (Relative to Remaining Service Area Mean)",
  #        y = "Cumulative density",
  #        title = "Pooled ECDFs of Remaining and Lost Wetlands in Mitigation Bank Service Areas\n(weighted by Lost Wetland Area)") +
  #   scale_color_manual(values = c(lost = "red", remaining = "blue"),
  #                      labels = c(lost = "Wetland Loss (1985–2021)",
  #                                 remaining = "Remaining Wetland (2021)"),
  #                      name = "") +
  #   theme_bw(base_size = 14) +
  #   theme(
  #     legend.position = "bottom",
  #     legend.title = element_text(size = 14),
  #     legend.text = element_text(size = 12),
  #     axis.title = element_text(size = 16),
  #     axis.text = element_text(size = 14),
  #     plot.title = element_text(size = 16),
  #     panel.grid.major = element_blank(),
  #     panel.grid.minor = element_blank()
  #   )+
  #   geom_vline(xintercept = 1, linetype = "dashed", colour = "gray40", linewidth = .8)
  
  # combined_plot <- pdf_plot / ecdf_plot
  
  message("All done! Returning plots and data.")
  list(
    pdf_plot      = pdf_plot,
    #ecdf_plot     = ecdf_plot,
    #combined_plot = combined_plot,
    data          = results_ecdf_full,
    p_value       = p_val,
    sig_stars     = sig_stars
  )
}

bank_summary <- readRDS(pth("Extractions and Summaries", "bank_summary.rds"))
sa_summary   <- readRDS(pth("Extractions and Summaries", "sa_summary.rds"))
loss_summary <- readRDS(pth("Extractions and Summaries", "loss_summary.rds"))

# FP_L Name update
sa_summary$sa_id <- gsub("FP_L_", "FP_AndL_", sa_summary$sa_id)
loss_summary$sa_id <- gsub("FP_L_", "FP_AndL_", loss_summary$sa_id)

# Add Source column
bank_summary <- bank_summary %>%
  rename(sa_id = bank_name) %>%
  mutate(Source = "Bank")

sa_summary <- sa_summary %>%
  mutate(Source = "Service Area")

sa_summary_subset <- sa_summary %>%
  filter(sa_id %in% bank_summary$sa_id)

loss_summary_subset = loss_summary %>% 
  filter(sa_id %in% bank_summary$sa_id )
# Merge
combined_summary <- bind_rows(bank_summary, sa_summary_subset, loss_summary_subset)
saveRDS(combined_summary, pth("Extractions and Summaries", "full_summary.rds"))

test = readRDS(pth("Extractions and Summaries", "full_summary.rds"))
names = unique(test$sa_id)
# names = names[1:6]
res <- compare_remaining_vs_lost(datadir, bank_filter = names)
# res$combined_plot      # shows PDF + ECDF with bracket & sig-stars
res$pdf_plot




##############################################################################
######                                                                ########
######          FIGURE 2 (updated, pooled service areas)              ########
######                                                                ########
##############################################################################


compare_remaining_vs_lost_pooled <- function(datadir, bank_pool_map, bank_summary_path, sample_size = 10000, seed = 123, filter = NULL) {
  suppressPackageStartupMessages({
    library(tidyverse)
    library(scales)
    library(vctrs)
  })
  
  message("Setting up pooled analysis with bank_pool_map...")
  
  # Load bank summary to get valid bank names
  bank_summary <- readRDS(bank_summary_path)
  valid_sa_names <- unique(bank_summary$bank_name)  
  message(sprintf("Loaded %d valid bank names from bank_summary", length(valid_sa_names)))
  
  # Create mapping of pools to representative bank names (for loss files)
  pool_representatives <- bank_pool_map %>%
    group_by(pool_id) %>%
    slice(1) %>%  # Take first bank as representative
    ungroup()
  
  # Get list of banks that are in pools
  pooled_banks <- unique(bank_pool_map$bank_name)
  
  # Get all available individual SA files
  sa_files <- list.files(file.path(datadir, "SA ECDF Functions"), full.names = FALSE)
  sa_names <- sub("^ecdf_fn_(.*)\\.rds$", "\\1", sa_files)
  message(sprintf("Found %d total SA files", length(sa_names)))
  
  # Filter to valid SA names from bank_summary 
  sa_names <- intersect(sa_names, valid_sa_names)
  message(sprintf("Filtered to valid SA names: %d", length(sa_names)))
  
  # Remove pooled banks from individual processing
  individual_banks <- setdiff(sa_names, pooled_banks)
  
  if (!is.null(filter)) {
    individual_banks <- individual_banks[individual_banks %in% filter]
    # Also filter pools if any of their banks are in the filter
    relevant_pools <- bank_pool_map %>%
      filter(bank_name %in% filter) %>%
      pull(pool_id) %>%
      unique()
    pool_representatives <- pool_representatives %>%
      filter(pool_id %in% relevant_pools)
  }
  
  # Debug: Show what we're working with
  message(sprintf("Total SA files found: %d", length(sa_names)))
  message(sprintf("Banks in pools (to be removed): %d", length(pooled_banks)))
  message(sprintf("Individual banks after removing pooled: %d", length(individual_banks)))
  message(sprintf("Number of pools: %d", nrow(pool_representatives)))
  
  # Calculate total markets: individual banks + pools
  n_individual <- length(individual_banks)
  n_pools <- nrow(pool_representatives)
  total_markets <- n_individual + n_pools
  
  message(sprintf("Expected calculation: %d individual + %d pools = %d total markets", 
                  n_individual, n_pools, total_markets))
  
  # Check if we have the expected 936 markets
  if (total_markets != 936) {
    warning(sprintf("Expected 936 markets but found %d. Please verify:", total_markets))
    warning(sprintf("  - Original SA files: %d", length(sa_names)))
    warning(sprintf("  - Pooled banks to remove: %d", length(pooled_banks)))
    warning(sprintf("  - Should be: %d - %d + %d = 936", length(sa_names), length(pooled_banks), n_pools))
  }
  
  results_list <- vector("list", total_markets * 2)
  j <- 1
  successful_markets <- 0
  set.seed(seed)
  
  # Process individual banks (non-pooled)
  message("Processing individual (non-pooled) banks...")
  for (i in seq_along(individual_banks)) {
    bank <- individual_banks[i]
    market_num <- i
    message(sprintf("[%d/%d] Processing individual bank: %s", market_num, total_markets, bank))
    
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
    
    row_weight <- loss_n_obs / length(loss_data)
    
    remaining <- tibble(
      market_id = bank, 
      status = "remaining", 
      value = sa_data,
      normalized_value = sa_data / bank_remaining_mean,
      weight = row_weight
    ) %>% filter(is.finite(normalized_value) & normalized_value > 0)
    
    lost <- tibble(
      market_id = bank, 
      status = "lost", 
      value = loss_data,
      normalized_value = loss_data / bank_remaining_mean,
      weight = row_weight
    ) %>% filter(is.finite(normalized_value) & normalized_value > 0)
    
    results_list[[j]] <- remaining; j <- j + 1
    results_list[[j]] <- lost;      j <- j + 1
    successful_markets <- successful_markets + 1
  }
  
  # Process pooled banks
  message("Processing pooled banks...")
  for (i in seq_len(nrow(pool_representatives))) {
    pool_id <- pool_representatives$pool_id[i]
    representative_bank <- pool_representatives$bank_name[i]
    market_num <- n_individual + i
    
    message(sprintf("[%d/%d] Processing pool: %s (using %s for loss)", 
                    market_num, total_markets, pool_id, representative_bank))
    
    # Use pooled SA ECDF
    pooled_sa_file <- file.path(datadir, "Pooled HUC8 ECDFs", sprintf("ecdf_pool_vs_sa_%s.rds", pool_id))
    
    # Use representative bank for loss ECDF
    loss_file <- file.path(datadir, "SA ECDF Functions - Loss", sprintf("ecdf_fn_%s.rds", representative_bank))
    
    if (!file.exists(pooled_sa_file) || !file.exists(loss_file)) {
      message("  - Skipping – missing ECDF file(s)")
      next
    }
    
    # For pooled files, we want $sa_rem_ecdf$housing_units
    pooled_data <- readRDS(pooled_sa_file)
    sa_ecdf <- pooled_data$sa_rem_ecdf$housing_units
    
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
    
    pool_remaining_mean <- mean(sa_data)
    
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
    
    row_weight <- loss_n_obs / length(loss_data)
    
    remaining <- tibble(
      market_id = pool_id, 
      status = "remaining", 
      value = sa_data,
      normalized_value = sa_data / pool_remaining_mean,
      weight = row_weight
    ) %>% filter(is.finite(normalized_value) & normalized_value > 0)
    
    lost <- tibble(
      market_id = pool_id, 
      status = "lost", 
      value = loss_data,
      normalized_value = loss_data / pool_remaining_mean,
      weight = row_weight
    ) %>% filter(is.finite(normalized_value) & normalized_value > 0)
    
    results_list[[j]] <- remaining; j <- j + 1
    results_list[[j]] <- lost;      j <- j + 1
    successful_markets <- successful_markets + 1
  }
  
  if (successful_markets == 0) stop("No markets with both remaining and lost ECDFs were processed.")
  message("Successfully processed ", successful_markets, " market(s) out of ", total_markets, " expected")
  
  results <- bind_rows(results_list)
  rm(results_list); gc()
  
  message("Computing ECDF data …")
  results_ecdf_full <- results %>%
    group_by(market_id, status) %>%
    arrange(normalized_value) %>%
    mutate(ecdf_y = row_number() / n()) %>%
    ungroup()
  
  # Statistics for the PDF panel
  loss_mean <- weighted.mean(results$normalized_value[results$status == "lost"],
                             w = results$weight[results$status == "lost"], na.rm = TRUE)
  
  # t-test (unweighted)
  p_val <- t.test(normalized_value ~ status, data = results)$p.value
  
  # Create the plot using the separate function
  pdf_plot <- create_wetland_pdf_plot(results_ecdf_full, loss_mean, p_val)
  
  message("All done! Returning plots and data.")
  list(
    pdf_plot      = pdf_plot,
    data          = results_ecdf_full,
    p_value       = p_val,
    loss_mean     = loss_mean,
    n_markets     = successful_markets,
    n_individual  = n_individual,
    n_pools       = n_pools
  )
}


datadir <- "L:/Wetland Flood Mitigation/ECDF Functions"
bank_summary_path <- "CONUS/Extractions and Summaries/bank_summary.rds"
res <- compare_remaining_vs_lost_pooled(datadir, bank_pool_map, bank_summary_path, sample_size = 10000, seed = 123)
res$pdf_plot

# Check that we have 936 markets:
message(sprintf("Total markets processed: %d (Individual: %d, Pools: %d)",
                res$n_markets, res$n_individual, res$n_pools))


# CREATE PDF PLOT HERE:

create_wetland_pdf_plot <- function(results_data, loss_mean, p_val,
                                    bracket_pad   = 0.10,   # ← raise bracket
                                    star_pad_frac = 0.05) { # ← raise stars
  sig_stars <- case_when(
    p_val < 0.001 ~ "***",
    p_val < 0.01  ~ "**",
    p_val < 0.05  ~ "*",
    TRUE          ~ "ns"
  )
  
  dens_rem  <- density(results_data$normalized_value[results_data$status == "remaining"])
  dens_lost <- density(results_data$normalized_value[results_data$status == "lost"])
  y_max     <- max(dens_rem$y, dens_lost$y, na.rm = TRUE)
  
  y_bracket <- y_max * (1 + bracket_pad)           # was 1.05
  y_tick    <- y_max * (1 + bracket_pad - 0.03)    # keep ticks a bit lower
  
  x_left  <- min(1, loss_mean)
  x_right <- max(1, loss_mean)
  
  ggplot(results_data,
         aes(x = normalized_value,
             colour = status,
             weight = weight)) +
    stat_density(geom = "line", position = "identity", linewidth = 1.2) +
    geom_vline(xintercept = 1,         linetype = "dashed",
               colour = "blue",  linewidth = .8) +
    geom_vline(xintercept = loss_mean, linetype = "dashed",
               colour = "red",   linewidth = .8) +
    annotate("segment", x = x_left,  xend = x_right,
             y = y_bracket, yend = y_bracket) +
    annotate("segment", x = x_left,  xend = x_left,
             y = y_bracket, yend = y_tick) +
    annotate("segment", x = x_right, xend = x_right,
             y = y_bracket, yend = y_tick) +
    annotate("text",
             x = (x_left + x_right) / 2 - 0.2,
             y = y_bracket + star_pad_frac * y_max,
             label = sig_stars, vjust = 0) +
    scale_x_log10(
      limits = c(0.01, 100),
      breaks = sort(unique(c(0.1, 1, 10, 100, signif(loss_mean, 2)))),
      labels = function(x) format(x, digits = 2, trim = TRUE),
      expand  = c(0, 0)
    ) +
    expand_limits(y = y_bracket + (star_pad_frac + 0.01) * y_max) +
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
      legend.text   = element_text(size = 12),
      axis.title    = element_text(size = 16),
      axis.text     = element_text(size = 14),
      plot.title    = element_text(size = 16),
      panel.grid.major = element_blank(),
      panel.grid.minor = element_blank()
    )
}

higher_plot <- create_wetland_pdf_plot(res$data,
                                       loss_mean = res$loss_mean,
                                       p_val     = res$p_value,
                                       bracket_pad   = 0.25,   
                                       star_pad_frac = 0.12)   
print(higher_plot)

ggsave(file.path("CONUS/Graphics/pdfplot_pooled.png"),  plot_result$plot, width=12, height=8, dpi=300)

# -------------

##############################################################################
######                                                                ########
######          FIGURE 3 (old version, no pooled service areas)       ########
######                                                                ########
##############################################################################

suppressPackageStartupMessages({
  library(future)
  library(future.apply)
  library(tidyverse)
  library(progressr)
  library(data.table)
})

# ------------------------------------------------------------------
#                      Paths & helpers
# ------------------------------------------------------------------
datadir <- "C:/Users/indumati/Box/Paper2_final"      # project root
pth     <- function(...) file.path(datadir, ...)     # quick path builder

ecdf_dir   <- pth("Extractions and Summaries", "ECDF Functions")                  # top of ECDF tree
bank_dir   <- pth("Extractions and Summaries")  # where *_summary.rds live


# ------------------------------------------------------------------
#                       ECDF file lists
# ------------------------------------------------------------------
name_files <- list.files(pth(ecdf_dir, "SA ECDF Functions"), full.names = FALSE)
sa_names <- sub("^ecdf_fn_(.*)\\.rds$", "\\1", name_files)

bankfiles <- list.files(pth(ecdf_dir, "Bank ECDF Functions"),
                        pattern = "^bank_ecdf_fn_.*\\.rds$", full.names = TRUE)
safiles   <- list.files(pth(ecdf_dir, "SA ECDF Functions"),
                        pattern = "^ecdf_fn_.*\\.rds$", full.names = TRUE)
lossfiles <- list.files(pth(ecdf_dir, "Loss ECDF Functions"),
                        pattern = "^ecdf_fn_.*\\.rds$", full.names = TRUE)

# ------------------------------------------------------------------
#                       helper functions
# ------------------------------------------------------------------

get_ecdf_data <- function(ecdf_obj) {
  dat <- environment(ecdf_obj)$x
  if (length(dat) > 1e6) dat <- sample(dat, 1e6)
  dat
}

signed_area_diff <- function(ecdf2, ecdf1, grid_size = 1000) {
  x1 <- get_ecdf_data(ecdf1)
  x2 <- get_ecdf_data(ecdf2)
  grid <- seq(min(c(x1, x2)), max(c(x1, x2)), length.out = grid_size)
  delta_x <- diff(range(grid)) / (grid_size - 1)
  sum(ecdf1(grid) - ecdf2(grid)) * delta_x
}

quantile_diff <- function(x1, x2, probs = c(.25, .5, .75, .95, .99)) {
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

process_unit <- function(unit_name) {
  bankfile <- pth(ecdf_dir, "Bank ECDF Functions",        paste0("bank_ecdf_fn_", unit_name, ".rds"))
  safile   <- pth(ecdf_dir, "SA ECDF Functions",          paste0("ecdf_fn_",    unit_name, ".rds"))
  lossfile <- pth(ecdf_dir, "SA ECDF Functions - Loss",   paste0("ecdf_fn_",    unit_name, ".rds"))
  
  bank <- if (file.exists(bankfile)) readRDS(bankfile)$housing_units else NA
  sa   <- if (file.exists(safile))   readRDS(safile)$housing_units   else NA
  loss <- if (file.exists(lossfile)) readRDS(lossfile)$housing_units else NA
  
  bind_rows(
    make_row(unit_name, "loss_vs_sa",   loss, sa),
    make_row(unit_name, "loss_vs_bank", loss, bank),
    make_row(unit_name, "sa_vs_bank",   sa,   bank)
  )
}

# ------------------------------------------------------------------
#                        run the comparisons
# ------------------------------------------------------------------
plan(multisession)
handlers(global = TRUE)

with_progress({
  p <- progressor(along = bank_names)
  results_df_long <- future_lapply(bank_names, function(n) { p(); process_unit(n) })
})
results_df_long <- bind_rows(results_df_long)
plan(sequential)

# ------------------------------------------------------------------
#         summaries & additional ECDF normalisation (unchanged)
# ------------------------------------------------------------------
bank_summary  <- readRDS(pth(bank_dir, "bank_summary.rds"))
sa_summary    <- readRDS(pth(bank_dir, "sa_summary.rds"))
loss_summary  <- readRDS(pth(bank_dir, "loss_summary.rds"))

sa_summary$sa_id   <- gsub("FP_L_", "FP_AndL_", sa_summary$sa_id)
loss_summary$sa_id <- gsub("FP_L_", "FP_AndL_", loss_summary$sa_id)




library(future)
library(future.apply)
library(ggplot2)
library(tidyverse)
library(tibble)
library(progressr)
library(data.table)
library(dplyr)

datadir="C:/Users/indumati/Box/Paper2_final/ECDF Functions"

#get unique bank names
names=list.files(paste0(datadir,"\\SA ECDF Functions"))
names = sub("^ecdf_fn_(.*)\\.rds$", "\\1", names)

#calculate differences in means, medians and area between ECDFs for:
#1. lost wetlands vs service areas
#2. lost wetlands vs banks
#3. service areas vs banks


# Functions
get_ecdf_data <- function(ecdf_obj) {
  dat=environment(ecdf_obj)$x
  if(length(dat)>1e6) dat=sample(dat,1e6,replace=FALSE)
  return(dat)
}

signed_area_diff <- function(ecdf2, ecdf1, grid_size = 1000) {
  x1 <- get_ecdf_data(ecdf1)
  x2 <- get_ecdf_data(ecdf2)
  x_range <- range(c(x1, x2))
  grid <- seq(x_range[1], x_range[2], length.out = grid_size)
  delta_x <- diff(range(grid)) / (grid_size - 1)
  sum(ecdf1(grid) - ecdf2(grid)) * delta_x
}

quantile_diff <- function(x1, x2, probs = c(0.25, 0.5, 0.75,0.95,0.99)) {
  q1 <- quantile(x1, probs = probs)
  q2 <- quantile(x2, probs = probs)
  setNames(q1 - q2, paste0("qdiff_", probs*100))
}

make_row <- function(unit_name, label, ecdf1, ecdf2) {
  if (!is.function(ecdf1) || !is.function(ecdf2)) {
    return(data.frame(
      unit = unit_name,
      comparison = label,
      mean_diff = NA,
      area_diff = NA,
      qdiff_25 = NA,
      qdiff_50 = NA,
      qdiff_75 = NA,
      qdiff_95 = NA,
      qdiff_99 = NA
    ))
  }
  
  x1 <- get_ecdf_data(ecdf1)
  x2 <- get_ecdf_data(ecdf2)
  qdiffs <- quantile_diff(x1, x2)
  
  data.frame(
    unit = unit_name,
    comparison = label,
    mean_diff = mean(x1) - mean(x2),
    area_diff = signed_area_diff(ecdf1, ecdf2),
    qdiff_25 = qdiffs["qdiff_25"],
    qdiff_50 = qdiffs["qdiff_50"],
    qdiff_75 = qdiffs["qdiff_75"],
    qdiff_95 = qdiffs["qdiff_95"],
    qdiff_99 = qdiffs["qdiff_99"]
  )
}

bankfiles=list.files(paste0(datadir,"Extractions and Summaries/Bank ECDF Functions"),full.names=TRUE)
safiles=list.files(paste0(datadir,"Extractions and Summaries/SA ECDF Functions"),full.names=TRUE)
lossfiles=list.files(paste0(datadir,"Extractions and Summaries/Loss ECDF Functions"),full.names=TRUE)


x_grid <- 10^seq(-2, 3.5, length.out = 750)
#next compare: mean bank values (i.e. reasonable mitigation) against different quantiles of loss distribution
#look at ratios - i.e. what is ratio of flood benefit at different parts of the wetland loss distribution?
process_unit <- function(unit_name, datadir=datadir, bankfiles=bankfiles, lossfiles=lossfiles, x_grid=x_grid) {
  # File paths
  bankfile <- paste0(datadir, "Extractions and Summaries/Bank ECDF Functions/bank_ecdf_fn_", unit_name, ".rds")
  lossfile <- paste0(datadir, "Extractions and Summaries/Loss ECDF Functions/ecdf_fn_", unit_name, ".rds")
  
  # Load data if files exist
  bank_data <- if (bankfile %in% bankfiles) readRDS(bankfile) else NA
  loss_data <- if (lossfile %in% lossfiles) readRDS(lossfile) else NA
  
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


plan(multisession)  
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

fwrite(results_long_normalized,file=paste0(datadir,"\\normalized_ecdfs_loss.csv"))

# Highest number of loss cells:
loss_cell_counts <- results_long_normalized %>%
  filter(comparison == "normalized_loss_units") %>%
  group_by(unit) %>%
  slice(1) %>%             # take the first row per unit
  ungroup() %>%
  dplyr::select(unit, n_obs) %>%
  arrange(desc(n_obs))
loss_cell_counts

######## cutoff thresholds  ------------------

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


# ------------------------

# 1. Define the reference x-values
ref_values <- c(1, 1.5, 2)

# 2. Interpolate the mean CDF at each reference for each comparison
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

# 3. Add horizontal lines & labels to your existing plot ‘a’
#    We anchor the left side at x = 0.1 for readability in log-scale plots:

# Then update your plot:
a +
  geom_vline(
    xintercept = 1,
    linetype = "dashed",
    color = "gray20"
  ) +
  geom_vline(
    xintercept = 1.5,
    linetype = "dashed",
    color = "gray50"
  ) +
  geom_vline(
    xintercept = 2,
    linetype = "dashed",
    color = "gray80"
  ) +
  geom_segment(
    data = crossings,
    aes(x = xmin_val, xend = ref_x, y = cdf_value, yend = cdf_value),
    inherit.aes = FALSE,
    linetype = "dashed",
    color = "black"
  ) +
  scale_x_log10(labels = scales::label_number())


# Filter only housing unit rows
results_hu <- results_long_normalized %>%
  filter(comparison == "normalized_loss_units")

# Recalculate ecdf_summary for housing units only
ecdf_summary_hu <- results_hu %>%
  group_by(x, comparison) %>%
  summarize(mean_y = weighted.mean(y, w = n_obs), .groups = "drop")

# big_cypress_line <- results_hu %>% filter(unit == "Big_Cypress_MB_Phase_I-V")
# 
# big_cypress_label_point <- big_cypress_line %>%
#   filter(y >= 0.4 & y <= 0.6) %>%
#   slice_min(x)

# Assign colors to the reference lines
ref_line_colors <- c("gray10", "gray25", "gray40")
names(ref_line_colors) <- c("below_1", "below_1_5", "below_2")

# Merge colors into crossings
crossings <- crossings %>%
  mutate(line_color = ref_line_colors[key])
crossings_hu <- crossings %>% filter(comparison == "normalized_loss_units")

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


########################################################
########                                          ######
########  Figure 3 (updated version, pooled SAs)  ######
########                                          ######
########################################################


# ──────────────────────────────────────────────────────────────
# 3.  Helper functions  (unchanged from your script)
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

# Pool representative workflow with exact same process_unit logic
library(dplyr)
library(data.table)
library(future)
library(future.apply)
library(progressr)
datadir = base_dir <- "C:\\Users\\indumati\\Box\\ECDF Functions"
# Setup directories and file lists
base_dir <- "C:\\Users\\indumati\\Box\\ECDF Functions"
bankfiles <- list.files(file.path(base_dir, "Bank ECDF Functions"), full.names = TRUE)
lossfiles <- list.files(file.path(base_dir, "SA ECDF Functions - Loss"), full.names = TRUE)

# Get available unit names from bank files
bank_names <- list.files(file.path(base_dir, "Bank ECDF Functions"))
bank_names <- sub("^bank_ecdf_fn_(.*)\\.rds$", "\\1", bank_names)

loss_names <- list.files(file.path(base_dir, "SA ECDF Functions - Loss"))
loss_names <- sub("^ecdf_fn_(.*)\\.rds$", "\\1", loss_names)

# Units that have both bank AND loss files 
complete_units <- intersect(bank_names, loss_names)
cat("Units with both bank and loss files:", length(complete_units), "\n")

bank_pool_map = readRDS("C:\\Users\\indumati\\Box\\Paper2_final\\Extractions and Summaries\\ECDF Functions\\bank_pool_map.rds")
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

datadir = "C:\\Users\\indumati\\Box\\Paper2_final"
# process_unit <- function(unit_name, datadir=datadir, bankfiles=bankfiles, lossfiles=lossfiles, x_grid=x_grid) {
#   # Always load loss data
#   lossfile <- paste0(datadir, "ECDF Functions/SA ECDF Functions - Loss/ecdf_fn_", unit_name, ".rds")
#   loss_data <- if (lossfile %in% lossfiles) readRDS(lossfile) else NA
# 
#   # Determine if unit is in a pool
#   pool_id <- pool_lookup$pool_id[pool_lookup$sa_id == unit_name]
# 
#   if (length(pool_id) > 0) {
#     # POOLED UNIT: Load from pooled directory
#     poolfile <- paste0(datadir, "ECDF Functions/Pooled HUC8 ECDFs/ecdf_pool_vs_sa_", pool_id[1], ".rds")
#     if (file.exists(poolfile)) {
#       pool_data <- readRDS(poolfile)
#       bank_data <- pool_data$bank_ecdf
#     } else {
#       return(NULL)
#     }
#   } else {
#     # SOLO UNIT: Load from individual bank directory
#     bankfile <- paste0(datadir, "ECDF Functions/Bank ECDF Functions/bank_ecdf_fn_", unit_name, ".rds")
#     bank_data <- if (bankfile %in% bankfiles) readRDS(bankfile) else NA
#   }
# 
#   # Validate required objects
#   if (!is.list(bank_data) || !is.list(loss_data)) return(NULL)
#   if (!is.function(bank_data$housing_units) || !is.function(loss_data$housing_units)) return(NULL)
#   if (!is.function(bank_data$housing_value) || !is.function(loss_data$housing_value)) return(NULL)
# 
#   # Extract raw values
#   bank_vals  <- environment(bank_data$housing_units)$x
#   bank_hval  <- environment(bank_data$housing_value)$x
#   loss_vals  <- environment(loss_data$housing_units)$x
#   loss_hval  <- environment(loss_data$housing_value)$x
# 
#   # Compute means
#   bank_mean_vals <- mean(bank_vals, na.rm = TRUE)
#   bank_mean_hval <- mean(bank_hval, na.rm = TRUE)
# 
#   # Normalize accordingly
#   loss_vals_norm <- loss_vals / bank_mean_vals
#   loss_hval_norm <- loss_hval / bank_mean_hval
# 
#   n_obs <- length(loss_vals_norm)
# 
#   # Optional downsampling
#   if (n_obs > 1e5) {
#     sampled_indices <- sample(seq_len(n_obs), 1e5)
#     loss_vals_norm <- loss_vals_norm[sampled_indices]
#     loss_hval_norm <- loss_hval_norm[sampled_indices]
#   }
# 
#   # Evaluate ECDFs
#   ecdf_vals  <- ecdf(loss_vals_norm)
#   ecdf_hvals <- ecdf(loss_hval_norm)
# 
#   # Output both sets as tibbles
#   result_units <- tibble(
#     unit = unit_name,
#     comparison = "normalized_loss_units",
#     x = x_grid,
#     y = ecdf_vals(x_grid),
#     n_obs = n_obs
#   )
# 
#   result_values <- tibble(
#     unit = unit_name,
#     comparison = "normalized_loss_values",
#     x = x_grid,
#     y = ecdf_hvals(x_grid),
#     n_obs = n_obs
#   )
# 
#   # Return combined result
#   bind_rows(result_units, result_values)
# }

process_unit <- function(unit_name, datadir=datadir, bankfiles=bankfiles, lossfiles=lossfiles, x_grid=x_grid) {
  # Always load loss data
  lossfile <- paste0(datadir, "ECDF Functions/SA ECDF Functions - Loss/ecdf_fn_", unit_name, ".rds")
  loss_data <- if (lossfile %in% lossfiles) readRDS(lossfile) else NA
  
  # Determine if unit is in a pool
  pool_id <- pool_lookup$pool_id[pool_lookup$sa_id == unit_name]
  
  if (length(pool_id) > 0) {
    # POOLED UNIT: Load from pooled directory
    poolfile <- paste0(datadir, "ECDF Functions/Pooled HUC8 ECDFs/ecdf_pool_vs_sa_", pool_id[1], ".rds")
    if (file.exists(poolfile)) {
      pool_data <- readRDS(poolfile)
      bank_data <- pool_data$bank_ecdf
    } else {
      return(NULL)
    }
  } else {
    # SOLO UNIT: Load from individual bank directory
    bankfile <- paste0(datadir, "ECDF Functions/Bank ECDF Functions/bank_ecdf_fn_", unit_name, ".rds")
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
plan(multisession, workers= 5)
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


fwrite(results_long_normalized,file=paste0(base_dir,"\\poolednormalized_ecdfs_loss.csv"))

results_long_normalized = 
results_hu <- results_long_normalized %>%
  filter(comparison == "normalized_loss_units")

# -----------------------------------------------------------------------------------------
# Get weighted mean -----------------------------------------------------------------------
# -----------------------------------------------------------------------------------------

# Use only the housing-units-normalized ECDFs
results_hu <- results_long %>%
  filter(comparison == "normalized_loss_units")

# 1) Get mean ratio per unit from its ECDF, then 2) take weighted mean across units
unit_means <- results_hu %>%
  arrange(unit, x) %>%
  group_by(unit) %>%
  mutate(delta = y - lag(y, default = 0)) %>%          # ECDF step sizes
  summarize(
    mean_ratio_unit = sum(x * delta, na.rm = TRUE),    # E[X | unit]
    n_obs = first(n_obs),                              # weight by observations in that unit
    .groups = "drop"
  )

mean_ratio_overall <- with(unit_means,
                           weighted.mean(mean_ratio_unit, w = n_obs, na.rm = TRUE)
)

mean_ratio_overall #4.148653

#------------------------------------------------------------

# Recalculate ecdf_summary for housing units only

ecdf_summary_hu <- results_hu %>%
  arrange(unit, x) %>%
  group_by(unit) %>%
  # ensure each unit has y defined for every x in the grid
  complete(x = x_grid) %>%
  arrange(x, .by_group = TRUE) %>%
  # carry forward step values; before first step, CDF = 0
  fill(y, .direction = "down") %>%
  mutate(y = if_else(is.na(y), 0, y)) %>%
  ungroup() %>%
  group_by(x) %>%
  summarize(
    mean_y = weighted.mean(y, w = n_obs, na.rm = TRUE),
    .groups = "drop"
  )
# ecdf_summary <- results_hu %>%
#   group_by(x, comparison) %>%
#   summarize(mean_y = weighted.mean(y, w = n_obs), .groups = "drop")
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

ggsave(file.path("Figures/ecdfplot_pooled.png"),  fig3_ecdf, width=12, height=8, dpi=300)

############################################
############ Figure 3 histogram ############
############################################

# ------------------ Libraries ------------------
library(tidyverse)   # dplyr, tibble, purrr, ggplot2, stringr, etc.

# ------------------ Paths ----------------------
# set this to your ECDF Functions folder (same as in your workflow)
datadir  <- "C:\\Users\\indumati\\Box\\Paper2_final"
bank_dir <- file.path(datadir, "Extractions and Summaries/ECDF Functions/Bank ECDF Functions")
loss_dir <- file.path(datadir, "Extractions and Summaries/ECDF Functions/SA ECDF Functions - Loss")

# ------------------ Helpers --------------------
# Extract raw x vector from an ecdf() function saved in RDS
get_ecdf_x <- function(fn) {
  if (!is.function(fn)) return(numeric(0))
  x <- environment(fn)$x
  x[is.finite(x)]
}

# Load a unit's RDS bundles; require HU, HV, and population to exist
load_unit <- function(unit) {
  bf <- file.path(bank_dir, paste0("bank_ecdf_fn_", unit, ".rds"))
  lf <- file.path(loss_dir, paste0("ecdf_fn_", unit, ".rds"))
  if (!file.exists(bf) || !file.exists(lf)) return(NULL)
  bank <- readRDS(bf)
  loss <- readRDS(lf)
  req  <- c("housing_units", "housing_value", "population")
  if (!all(req %in% names(bank)) || !all(req %in% names(loss))) return(NULL)
  list(unit = unit, bank = bank, loss = loss)
}

# ------------------ Discover units ------------------
bankfiles <- list.files(bank_dir, full.names = TRUE)
lossfiles <- list.files(loss_dir, full.names = TRUE)

units_bank <- basename(bankfiles) |>
  str_replace("^bank_ecdf_fn_", "") |>
  str_replace("\\.rds$", "")

units_loss <- basename(lossfiles) |>
  str_replace("^ecdf_fn_", "") |>
  str_replace("\\.rds$", "")

units <- intersect(units_bank, units_loss)

# ------------------ Load bundles ------------------
bundles <- lapply(units, load_unit)
bundles <- Filter(Negate(is.null), bundles)

# ------------------ Build unit_means ------------------
# (includes housing units, housing value, and population; we'll plot HU)
unit_means <- purrr::map_dfr(bundles, function(b) {
  bu <- get_ecdf_x(b$bank$housing_units)
  lu <- get_ecdf_x(b$loss$housing_units)
  bv <- get_ecdf_x(b$bank$housing_value)
  lv <- get_ecdf_x(b$loss$housing_value)
  bp <- get_ecdf_x(b$bank$population)
  lp <- get_ecdf_x(b$loss$population)
  
  tibble(
    unit = b$unit,
    # housing units
    mean_bank_units = mean(bu, na.rm = TRUE),
    mean_loss_units = mean(lu, na.rm = TRUE),
    n_bank_units    = length(bu),
    n_loss_units    = length(lu),
    # housing value
    mean_bank_value = mean(bv, na.rm = TRUE),
    mean_loss_value = mean(lv, na.rm = TRUE),
    n_bank_value    = length(bv),
    n_loss_value    = length(lv),
    # population
    mean_bank_pop   = mean(bp, na.rm = TRUE),
    mean_loss_pop   = mean(lp, na.rm = TRUE),
    n_bank_pop      = length(bp),
    n_loss_pop      = length(lp)
  )
})

# ------------------ Histogram (housing units) ------------------
library(ggplot2)
library(dplyr)

hist_df <- unit_means %>%
  filter(is.finite(mean_bank_units), mean_bank_units > 0) %>%
  mutate(ratio_loss_bank = mean_loss_units / mean_bank_units)

p_hist_hu <- ggplot(hist_df, aes(x = ratio_loss_bank)) +
  # histogram with finer bins
  geom_histogram(
    binwidth = 0.75,
    fill     = "#4C9BE8",
    color    = "white",
    alpha    = 0.8
  ) +
  # reference lines at 1 and 2
  geom_vline(xintercept = 1, linetype = "dashed", linewidth = 0.8, colour = "gray20") +
  # axis labels
  labs(
    x = "Flood Protection Index of Developed Wetlands, Relative to Bank Mean",
    y = "Number of Markets"
  ) +
  # custom breaks, clean integers
  scale_x_continuous(
    breaks = c(0, 1, 2, 10, 20, 50),
    expand = c(0, 0)   # no gap at left
  ) +
  scale_y_continuous(
    expand = c(0, 0)   # no gap at bottom
  ) +
  theme_bw(base_size = 14) +
  theme(
    axis.title       = element_text(size = 12),
    axis.text        = element_text(size = 12),
    axis.text.x      = element_text(angle = 45, hjust = 1),
    panel.grid.major = element_blank(),
    panel.grid.minor = element_blank(),
    # keep only axis lines
    panel.border     = element_blank(),
    axis.line        = element_line(color = "black")
  )

p_hist_hu
ggsave(file.path("CONUS/Graphics/fig3_hist.png"), p_hist_hu, width = 9, height = 5, dpi = 300)




# Pick ONE bin spec to use everywhere
bw <- 0.75
xmax_all <- max(hist_df$ratio_loss_bank, na.rm = TRUE)
# Explicit breaks so both plots bin identically
brks <- seq(0, xmax_all + bw, by = bw)

# --- MAIN ---
p_hist_hu <- ggplot(hist_df, aes(x = ratio_loss_bank)) +
  geom_histogram(
    breaks = brks,             # << identical breaks
    fill = "#4C9BE8", color = "white", alpha = 0.8
  ) +
  geom_vline(xintercept = 1, linetype = "dashed", linewidth = 0.8, colour = "gray20") +
  labs(x = "Flood Protection Index of Developed Wetlands, Relative to Bank Mean",
       y = "Number of Markets") +
  scale_x_continuous(breaks = c(0, 1, 2, 10, 20, 50), expand = c(0, 0)) +
  scale_y_continuous(expand = c(0, 0)) +
  theme_bw(base_size = 14) +
  theme(axis.text.x = element_text(angle = 45, hjust = 1),
        panel.grid.major = element_blank(),
        panel.grid.minor = element_blank(),
        panel.border = element_blank(),
        axis.line = element_line(color = "black"))

# --- INSET (zoom 0–20) ---
library(cowplot)

p_hist_inset <- ggplot(hist_df, aes(x = ratio_loss_bank)) +
  geom_histogram(
    breaks = brks,             # << same breaks as main
    fill = "#4C9BE8", color = "white", alpha = 0.9
  ) +
  geom_vline(xintercept = 1, linetype = "dashed", linewidth = 0.6, colour = "gray20") +
  labs(x = NULL, y = NULL) +
  coord_cartesian(xlim = c(0, 15)) +                    # zoom only, same bins
  scale_x_continuous(breaks = c(0, 1, 2, 4, 6, 10, 15, 20), expand = c(0, 0)) +
  scale_y_continuous(expand = expansion(mult = c(0, 0.05))) +  # small headroom
  theme_bw(base_size = 10) +
  theme(axis.text.x = element_text(angle = 45, hjust = 1),
        axis.text.y = element_text(size = 8),
        panel.grid.major = element_blank(),
        panel.grid.minor = element_blank(),
        panel.border = element_rect(color = "black", linewidth = 0.4, fill = NA),
        plot.margin = margin(2, 2, 2, 2))

# Overlay inset
p_hist_with_inset <- ggdraw(p_hist_hu) +
  draw_plot(p_hist_inset, x = 0.62, y = 0.55, width = 0.35, height = 0.35)

p_hist_with_inset

# ggsave("CONUS/Graphics/fig3_hist.png", p_hist_with_inset, width = 9, height = 5, dpi = 300)


## Hist with log scale

hist_df <- unit_means %>%
  filter(is.finite(mean_bank_units), mean_bank_units > 0, 
         is.finite(mean_loss_units), mean_loss_units > 0) %>%
  mutate(ratio_loss_bank = mean_loss_units / mean_bank_units)

cat("hist_df has", nrow(hist_df), "rows\n")

p_hist_hu <- ggplot(hist_df, aes(x = ratio_loss_bank)) +
  geom_histogram(
    bins = 50,  # Using bins instead of binwidth for log scale
    fill = "#4C9BE8",
    color = "white",
    alpha = 0.8
  ) +
  geom_vline(xintercept = 1, linetype = "dashed", linewidth = 0.8, colour = "gray20") +
  labs(
    x = "Flood Protection Index of Developed Wetlands, Relative to Bank Mean",
    y = "Number of Markets"
  ) +
  scale_x_log10(
    breaks = c(0.1, 1, 1.5, 2, 10, 100, 1000),
    labels = c("0.1", "1", "1.5", "2", "10", "100", "1000"),
    expand = c(0, 0)
  ) +
  scale_y_continuous(expand = c(0, 0)) +
  theme_bw(base_size = 14) +
  theme(
    axis.title = element_text(size = 12),
    axis.text = element_text(size = 12),
    axis.text.x = element_text(angle = 45, hjust = 1),
    panel.grid.major = element_blank(),
    panel.grid.minor = element_blank(),
    panel.border = element_blank(),
    axis.line = element_line(color = "black")
  )

p_hist_hu

################################# 
# SUMMARY STATS #################
#################################

### Averages 

# Define thresholds
cut_pts <- c(`gt1` = 1, `gte1.5` = 1.5, `gte2` = 2, `gte4.5` = 4.5, `gt10` = 10)

# Compute tail shares per unit directly from raw data
shares_per_unit <- map_dfr(names, function(nm) {
  loss_pth <- file.path(sa_loss_ecdf_dir, paste0("ecdf_fn_", nm, ".rds"))
  if (!file.exists(loss_pth)) return(NULL)
  
  bank_ecdf <- get_bank_ecdf(nm)
  if (is.null(bank_ecdf)) return(NULL)
  
  if (!is.function(bank_ecdf$housing_units)) return(NULL)
  loss_ecdf <- readRDS(loss_pth)
  if (!is.function(loss_ecdf$housing_units)) return(NULL)
  
  bank_vals <- environment(bank_ecdf$housing_units)$x
  loss_vals <- environment(loss_ecdf$housing_units)$x
  
  if (!is.numeric(bank_vals) || !is.numeric(loss_vals)) return(NULL)
  if (!length(bank_vals) || !length(loss_vals)) return(NULL)
  
  bank_mean <- mean(bank_vals, na.rm = TRUE)
  if (!is.finite(bank_mean) || bank_mean <= 0) return(NULL)
  
  # Normalize loss values by bank mean
  loss_ratio <- loss_vals / bank_mean
  
  # Calculate tail shares for each threshold
  shares <- map_dbl(cut_pts, function(thresh) {
    mean(loss_ratio > thresh, na.rm = TRUE) * 100
  })
  
  tibble(
    unit = nm,
    n_obs = length(loss_vals),
    !!!setNames(as.list(shares), names(cut_pts))
  )
})

# Overall weighted mean across all units
overall_wtd <- shares_per_unit %>%
  select(-unit) %>%
  summarise(across(-n_obs, ~ weighted.mean(.x, w = n_obs, na.rm = TRUE)))

print(overall_wtd)

# For a nice formatted table
percent_above_tbl <- tibble(
  threshold = names(cut_pts),
  x = cut_pts,
  pct_above = c(overall_wtd$gt1, overall_wtd$`gte1.5`, overall_wtd$gte2, 
                overall_wtd$`gte4.5`, overall_wtd$gt10)
) %>%
  mutate(pct_above = round(pct_above, 1))

print(percent_above_tbl)


###

######
# --- Inputs 
datadir          <- "L:/Wetland Flood Mitigation/ECDF Functions"
bank_ecdf_dir    <- file.path(datadir, "Bank ECDF Functions")
sa_loss_ecdf_dir <- file.path(datadir, "SA ECDF Functions - Loss")
pooled_dir       <- file.path(datadir, "Pooled HUC8 ECDFs")
bank_pool_map    <- readRDS(paste0(datadir, "/bank_pool_map.rds"))

# Pool lookup: pool_id per service-area id 
pool_lookup <- bank_pool_map %>% distinct(pool_id, sa_id = bank_name)

# Build unit list: one rep service area per pool + all solo service areas
bank_units  <- sub("^bank_ecdf_fn_(.*)\\.rds$", "\\1", list.files(bank_ecdf_dir))
loss_units  <- sub("^ecdf_fn_(.*)\\.rds$", "\\1", list.files(sa_loss_ecdf_dir))
complete    <- intersect(bank_units, loss_units)

# All service areas that are part of any pool
all_pooled_sa <- pool_lookup$sa_id

# One representative service area per pool (from those with complete data)
pools_with_complete <- pool_lookup %>%
  filter(sa_id %in% complete) %>%
  group_by(pool_id) %>% 
  slice(1) %>% 
  ungroup()
sa_multi <- pools_with_complete$sa_id

# Solo service areas: not in any pool
sa_solo <- setdiff(complete, all_pooled_sa)

# Final list of service area IDs to process
names <- c(sa_multi, sa_solo)

cat("Service areas representing pools:", length(sa_multi), "\n")
cat("Solo service areas:", length(sa_solo), "\n")
cat("Total service areas to process:", length(names), "\n")

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

# Area-weighted average across all service areas
area_weighted_avg <- with(ratio_df, weighted.mean(ratio_mean, w = n_pixels, na.rm = TRUE))
max_row <- ratio_df %>% arrange(desc(ratio_mean)) %>% slice(1)

p50_ratio <- median(ratio_df$ratio_mean, na.rm = TRUE)
p95_ratio <- quantile(ratio_df$ratio_mean, 0.95, na.rm = TRUE)

cat(glue(
  "wetlands lost to development exist near prior developed areas and therefore ",
  "provide relatively high levels of downstream flood protection compared to ",
  "wetlands created in compensation--on average ≈ {round(area_weighted_avg, 1)} times as much, ",
  "though in some cases up to {round(max_row$ratio_mean, 0)} times."
))

cat("\n\nDetails:\n")
cat("Number of service areas analyzed:", nrow(ratio_df), "\n")
cat("Median ratio:", round(p50_ratio, 1), "\n")
cat("95th percentile ratio:", round(p95_ratio, 1), "\n")
max_row


################# 
### Figure 4 ####
#################
# ECDF Lox and Cypress HU



# ECDFs --------------------------
ecdf_loss = readRDS(paste0(datadir,"ECDF Functions/SA ECDF Functions - Loss/ecdf_fn_Loxahatchee_MB.rds"))
ecdf_sa = readRDS(paste0(datadir, "ECDF Functions/SA ECDF Functions\\ecdf_fn_Loxahatchee_MB.rds"))
ecdf_bank = readRDS(paste0(datadir, "ECDF Functions/Bank ECDF Functions\\bank_ecdf_fn_Loxahatchee_MB.rds"))


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

# Plot: Only Loss vs Service Area if Bank is missing

df_hu = df_all %>% 
  filter(variable == "Housing Units")
ecdfplot_hu = ggplot(df_hu, aes(x = value, y = ecdf, color = source)) +
  geom_line(linewidth = 1.2) +
  scale_x_continuous(
    labels = comma,
    expand = c(0, 0)
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
    title = "Loxahatchee Mitigation Bank",
    x = "Downstream Housing Units",
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
    plot.margin = margin(t = 10, r = 50, b = 10, l = 10)
  )
ecdfplot_hu


# ECDF CYPRESS --------------------------------

ecdf_loss = readRDS("ECDF Functions/SA ECDF Functions - Loss/ecdf_fn_Big_Cypress_MB_Phase_I-V.rds")
ecdf_sa = readRDS("ECDF Functions/SA ECDF Functions/ecdf_fn_Big_Cypress_MB_Phase_I-V.rds")
ecdf_bank = readRDS("ECDF Functions/Bank ECDF Functions/bank_ecdf_fn_Big_Cypress_MB_Phase_I-V.rds")


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

# Plot: Only Loss vs Service Area if Bank is missing
df_hu = df_all %>% 
  filter(variable=="Housing Units")

library(scales)

ecdfplot_hu = ggplot(df_hu, aes(x = value, y = ecdf, color = source)) +
  geom_line(linewidth = 1.2) +
  scale_x_continuous(
    limits = c(0, 120000), 
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
    x = "Downstream Housing Units",
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

ecdfplot_hu


# Maps, Lox and Cypress HU --------------------

# LOX 
purrr::walk(c("future", "future.apply", "ggplot2", "tidyverse", "tibble", "progressr", 
              "data.table", "dplyr", "terra", "tigris", "sf", "scales", "ggspatial", 
              "rosm", "raster", "rnaturalearth"), library, character.only = TRUE)

setwd("L:\\Wetland Flood Mitigation")

loss_summary = readRDS("Extractions and Summaries/Extract Service Area - Loss/sa_summary_loss.rds")
bank_summary = readRDS("Extractions and Summaries/Extract Bank/bank_summary_2021.rds")

sas= readRDS("Service Areas/ServiceAreas_agg.rds")
banks = readRDS("Bank Footprints/footprints_and_buffers.rds")


#Find which states intersect this SA ----------------------------------------
# lox = test %>% 
# filter(sa_id=="Loxahatchee_MB")

lox_banks = banks[banks$Name =="Loxahatchee_MB"]
lox_sa = sas[sas$ID == "Loxahatchee_MB",]

lox_mosaic_loss = rast(paste0(datadir, "Service Area Mosaics/Loss Mosaics/Loxahatchee_MB_lossmosaic.tif"))
lox_mosaic_4326 = terra::project(lox_mosaic_loss, crs("EPSG:4326"))
states <- ne_states(country = "United States of America", returnclass = "sf")
states_albers <- st_transform(states, crs(lox_sa))
states_vect <- vect(states_albers)
intersected <- terra::intersect(states_vect, lox_sa)
intersected$name

fl_sf <- subset(states_albers, name %in% "Florida")

mosaic_df_loss <- as.data.frame(lox_mosaic_loss, xy = TRUE, na.rm = TRUE)
colnames(mosaic_df_loss)[3:5] <- c("housing_units", "housing_value", "population")

lox_sa <- st_as_sf(lox_sa)
lox_banks = st_as_sf(lox_banks)

# Note: ggspatial assumes your map is in lon/lat (EPSG:4326), so reproject everything
target_crs <- 4326

# Reproject everything to WGS84 (EPSG:4326)
lox_sa_4326 <- st_transform(lox_sa, target_crs)
lox_bank_4326 <- st_transform(lox_banks, target_crs)
fl_sf_4326 <- st_transform(fl_sf, target_crs)
coords_4326 <- terra::project(cbind(mosaic_df_loss$x, mosaic_df_loss$y), crs(lox_mosaic_loss), "EPSG:4326")

# Step 3: Create a new data frame with reprojected coordinates
mosaic_df_4326_loss <- data.frame(
  x = coords_4326[,1],
  y = coords_4326[,2],
  housing_units = mosaic_df_loss$housing_units,
  housing_value = mosaic_df_loss$housing_value,
  population = mosaic_df_loss$population
)

# Log transform for skew if needed
# mosaic_df_4326$log_housing_value <- log1p(mosaic_df_4326$housing_value)


# --------------------Before mapping, get bounding boxes -------------------------------

# 2) Get bounding box
bb <- st_bbox(lox_sa_4326)
# Optionally pad a bit if you want
bb["xmin"] <- bb["xmin"] - 0.2
bb["xmax"] <- bb["xmax"] + 0.2
bb["ymin"] <- bb["ymin"] - 0.2
bb["ymax"] <- bb["ymax"] + 0.2

# 3) Use xlim and ylim
xlim <- c(bb["xmin"], bb["xmax"])
ylim <- c(bb["ymin"], bb["ymax"])

# ------------------- Plot  -------------------------------

# --------- First, get bank average value in 2021 -----------------------------------------------
# Load the raster
bank_fill <- function(raster_path, polygon_sf, variable = "housing_units") {
  # Load raster and reproject to EPSG:4326
  raster_orig <- rast(raster_path)
  raster_4326 <- terra::project(raster_orig, "EPSG:4326")
  
  # Extract values with coordinates
  raster_df <- as.data.frame(raster_orig, xy = TRUE, na.rm = TRUE)
  colnames(raster_df)[3:5] <- c("housing_units", "housing_value", "population")
  
  # Reproject coordinates to EPSG:4326
  coords_4326 <- terra::project(cbind(raster_df$x, raster_df$y), crs(raster_orig), "EPSG:4326")
  
  # Combine into spatial points dataframe
  raster_df <- data.frame(
    x = coords_4326[,1],
    y = coords_4326[,2],
    housing_units = raster_df$housing_units,
    housing_value = raster_df$housing_value,
    population = raster_df$population
  )
  
  raster_points <- st_as_sf(raster_df, coords = c("x", "y"), crs = st_crs(polygon_sf))
  
  # Find points within polygon
  pts_in_poly <- raster_points[polygon_sf, , op = st_within]
  
  # Calculate mean of selected variable
  avg_val <- mean(pts_in_poly[[variable]], na.rm = TRUE)
  
  # Return the polygon with new variable added
  polygon_sf %>% mutate(!!variable := avg_val)
}

lox_bank_filled <- bank_fill(
  raster_path = paste0(datadir, "Service Area Mosaics/Loxahatchee_MB_mosaic.tif"),
  polygon_sf = lox_bank_4326,
  variable = "housing_units"  # can also be "housing_value" or "population"
)

lox_loss_filled = bank_fill(
  raster_path = paste0(datadir, "Service Area Mosaics/Loss Mosaics/Loxahatchee_MB_lossmosaic.tif"),
  polygon_sf = lox_sa_4326,
  variable = "housing_units"  # can also be "housing_value" or "population"
)

# ---------------------- Same LOSS map WITHOUT basemap ----------------------

cell_res <- res(lox_mosaic_4326)  # returns a vector: c(xres, yres)
lon_deg <- cell_res[1]
lat_deg <- cell_res[2]

places_fl <- places(state = "FL", year = 2021)  # or latest year available


lossplot_nbm = ggplot() +
  # Draw ocean *first* with soft blue
  annotate(
    "rect",
    xmin = xlim[1], xmax = xlim[2],
    ymin = ylim[1], ymax = ylim[2],
    fill = "lightblue3"
  ) +
  
  # Then draw land on top
  geom_sf(data = states, fill = "gray98", color = "gray70", linewidth = 0.2) +
  
  # Service area outline (legend label "Service Area")
  geom_sf(
    data  = lox_sa_4326,
    aes(color = "Service Area"),
    fill = NA,
    linewidth = 1
  ) +
  
  # Bank outline (legend label "Bank, avg. value 2021") + fill by housing_units
  # Fill only, no legend
  geom_sf(
    data = lox_bank_filled,
    aes(fill = housing_units),
    color = NA
  ) +
  
  # Outline only, for color legend
  geom_sf(
    data = lox_bank_filled,
    aes(color = "Bank, avg. value 2021"),
    fill = NA,
    linewidth = 0.4
  ) +
  
  # Downstream housing units (raster layer)
  geom_tile(
    data   = mosaic_df_4326_loss,  # Use the reprojected data
    aes(x = x, y = y, fill = housing_units),
    # Calculate appropriate width/height instead of using undefined variables
    width  = lon_deg,  # Example width in degrees (adjust as needed)
    height = lat_deg,  # Example height in degrees (adjust as needed)
    alpha  = 0.8
  )+
  
  geom_sf_label(
    data = places_fl %>% filter(NAME == "Cape Coral"),
    aes(label = NAME),
    stat = "sf_coordinates",
    nudge_y = 0.08,
    size = 3,
    label.size = 0.3,
    label.padding = unit(1, "mm")
  ) +
  
  # Naples and Fort Myers, no nudging
  geom_sf_label(
    data = places_fl %>% filter(NAME %in% c("Miami", "Hialea", "Fort Lauderdale", "Pembroke Pines")),
    aes(label = NAME),
    stat = "sf_coordinates",
    size = 3,
    label.size = 0.2,
    label.padding = unit(1, "mm")
  )+
  
  # Fill scale for raster
  scale_fill_viridis_c(
    option = "viridis",
    name = "Downstream Housing Units",
    labels = scales::comma_format(),
    direction = 1
  ) +
  
  # Color scale for outlines
  scale_color_manual(
    name = NULL,
    values = c(
      "Bank, avg. value 2021" = "red",
      "Service Area"          = "steelblue"
    ),
    guide = guide_legend(
      override.aes = list(
        fill = "transparent",  # 👈 This removes the unwanted gray fill
        shape = 22,            # optional: ensures a square key if needed
        size  = 3,             # optional: controls size of the square
        linewidth = 1
      )
    )
  ) +
  guides(
    color = guide_legend(order = 1),
    fill  = guide_colorbar(order = 2)
  ) +
  
  coord_sf(
    xlim = xlim,
    ylim = ylim,
    expand = FALSE, 
    crs = st_crs(target_crs)
  )+
  annotation_scale(location = "bl", width_hint = 0.25) +
  annotation_north_arrow(location = "tr", which_north = "true",
                         pad_x = unit(0.2, "in"), pad_y = unit(0.3, "in"),
                         style = north_arrow_fancy_orienteering) +
  
  theme_minimal(base_size = 13) +
  
  labs(
    title    = "Loxahatchee Mitigation Bank Service Area",
    subtitle = "Lost Wetlands (1985–2021)",
    fill     = "Downstream Housing Units", 
    x        = NULL,
    y        = NULL
  ) +
  
  theme(
    panel.grid = element_blank(),  
    legend.position   = "right",
    plot.title        = element_text(size = 13, hjust = 0.5),
    plot.subtitle     = element_text(size = 12, hjust = 0.5),
    legend.title      = element_text(size = 11),
    legend.text       = element_text(size = 11)
  )

lossplot_nbm


#  -------------------------------- Plot remaining SAs mosaic -----------------------------------------------------

# Load the raster
lox_mosaic_sa <- rast(paste0(datadir, "Service Area Mosaics/Loxahatchee_MB_mosaic.tif"))

# Extract values as a data frame (in native CRS, meters)
mosaic_df_sa <- as.data.frame(lox_mosaic_sa, xy = TRUE, na.rm = TRUE)
colnames(mosaic_df_sa)[3:5] <- c("housing_units", "housing_value", "population")

# Reproject raster coordinates to EPSG:4326
coords_4326 <- terra::project(cbind(mosaic_df_sa$x, mosaic_df_sa$y), crs(lox_mosaic_sa), "EPSG:4326")

# Create reprojected data frame for plotting
mosaic_df_sa <- data.frame(
  x = coords_4326[,1],
  y = coords_4326[,2],
  housing_units = mosaic_df_sa$housing_units,
  housing_value = mosaic_df_sa$housing_value,
  population = mosaic_df_sa$population
)


#  set extent from a relevant sf layer
zoom_ext_4326 <- st_bbox(lox_sa_4326)
bb <- zoom_ext_4326
bb["xmin"] <- bb["xmin"] - 0.2
bb["xmax"] <- bb["xmax"] + 0.2
bb["ymin"] <- bb["ymin"] - 0.2
bb["ymax"] <- bb["ymax"] + 0.2
xlim <- c(bb["xmin"], bb["xmax"])
ylim <- c(bb["ymin"], bb["ymax"])



########### PLOT SAME SA map without basemap #############

saplot_nbm = ggplot() +
  annotate("rect", xmin = xlim[1], xmax = xlim[2],
           ymin = ylim[1], ymax = ylim[2],
           fill = "lightblue3") +
  geom_sf(data = states, fill = "gray98", color = "gray70", linewidth = 0.2) +
  
  geom_sf(
    data  = lox_sa_4326,
    aes(color = "Service Area"),
    fill = NA,
    linewidth = 1
  ) +
  
  
  geom_tile(
    data = mosaic_df_sa |>
      dplyr::mutate(
        x = round(x, 5),
        y = round(y, 5)
      ),
    aes(x = x, y = y, fill = housing_units),
    width = lon_deg,
    height = lat_deg,
    color = NA,
    linewidth = 0
  )+
  geom_sf(
    data = lox_bank_filled,
    aes(
      fill  = housing_units,
      color = "Bank, avg. value 2021"
    ),
    linewidth = 0.4
  ) +
  
  geom_sf_label(
    data = places_fl %>% filter(NAME %in% c("Miami", "Hialea", "Fort Lauderdale", "Pembroke Pines")),
    aes(label = NAME),
    stat = "sf_coordinates",
    size = 3,
    label.size = 0.2,
    label.padding = unit(1, "mm")
  )+
  
  geom_sf_label(
    data = places_fl %>% filter(NAME %in% c("Naples", "Fort Myers")),
    aes(label = NAME),
    stat = "sf_coordinates",
    size = 3,
    label.size = 0.2,
    label.padding = unit(1, "mm")
  ) +
  
  scale_fill_viridis_c(
    option = "viridis",
    name = "Downstream Housing Units",
    labels = scales::comma_format(),
    direction = 1,
    na.value = "gray98"
  ) +
  
  scale_color_manual(
    name = NULL,
    values = c(
      "Bank, avg. value 2021" = "red",
      "Service Area"          = "steelblue"
    ),
    guide = guide_legend(
      override.aes = list(
        fill = NA,
        linewidth = 1
      )
    )
  ) +
  
  guides(
    color = guide_legend(order = 1),
    fill  = guide_colorbar(order = 2)
  ) +
  
  coord_sf(
    xlim = xlim, ylim = ylim,
    expand = FALSE,
    crs = st_crs(lox_sa_4326),
    default_crs = st_crs(lox_sa_4326)
  ) +
  annotation_scale(location = "bl", width_hint = 0.25) +
  annotation_north_arrow(location = "tr", which_north = "true",
                         pad_x = unit(0.2, "in"), pad_y = unit(0.3, "in"),
                         style = north_arrow_fancy_orienteering) +
  theme_minimal(base_size = 13) +
  labs(
    title    = "Loxahatchee Mitigation Bank Service Area",
    subtitle = "Remaining Wetlands (2021)",
    fill     = "Downstream Housing Units",
    x        = NULL,
    y        = NULL
  ) +
  theme(
    panel.grid       = element_blank(),
    legend.position  = "right",
    plot.title       = element_text(size = 13, hjust = 0.5),
    plot.subtitle    = element_text(size = 12, hjust = 0.5),
    legend.title     = element_text(size = 11),
    legend.text      = element_text(size = 11)
  )
saplot_nbm




# Map CYPRESS

test = readRDS(paste0(datadir,"Extractions and Summaries/Extract Service Area - Loss/sa_summary_loss.rds"))
test_bank = readRDS(paste0(datadir,"Extractions and Summaries/Extract Bank/bank_summary.rds"))

sas= readRDS(paste0(datadir,"Service Areas/ServiceAreas_agg.rds"))
banks = readRDS(paste0(datadir, "Footprints/footprints_and_buffers.rds"))


#Find which states intersect this SA ----------------------------------------
# cypress = test %>% 
# filter(sa_id=="Big_Cypress_MB_Phase_I-V")

cypress_banks = banks[banks$Name %in% c("Big_Cypress_MB_Phase_I-V", "Big_Cypress_MB_Phase_VI")]
cypress_banks = terra::union(cypress_banks)
cypress_sa = sas[sas$ID == "Big_Cypress_MB_Phase_I-V",]

cypress_mosaic_loss = rast(paste0(datadir, "Service Area Mosaics/Loss Mosaics/Big_Cypress_MB_Phase_I-V_lossmosaic.tif"))

states <- ne_states(country = "United States of America", returnclass = "sf")
states_albers <- st_transform(states, crs(cypress_sa))
states_vect <- vect(states_albers)
intersected <- terra::intersect(states_vect, cypress_sa)
intersected$name

fl_sf <- subset(states_albers, name %in% "Florida")

mosaic_df_loss <- as.data.frame(cypress_mosaic_loss, xy = TRUE, na.rm = TRUE)
colnames(mosaic_df_loss)[3:5] <- c("housing_units", "housing_value", "population")

cypress_sa <- st_as_sf(cypress_sa)
cypress_banks = st_as_sf(cypress_banks)

# Note: ggspatial assumes your map is in lon/lat (EPSG:4326), so reproject everything
target_crs <- 4326

# Reproject everything to WGS84 (EPSG:4326)
cypress_sa_4326 <- st_transform(cypress_sa, target_crs)
cypress_bank_4326 <- st_transform(cypress_banks, target_crs)
fl_sf_4326 <- st_transform(fl_sf, target_crs)
coords_4326 <- terra::project(cbind(mosaic_df_loss$x, mosaic_df_loss$y), crs(cypress_mosaic_loss), "EPSG:4326")

# Step 3: Create a new data frame with reprojected coordinates
mosaic_df_4326_loss <- data.frame(
  x = coords_4326[,1],
  y = coords_4326[,2],
  housing_units = mosaic_df_loss$housing_units,
  housing_value = mosaic_df_loss$housing_value,
  population = mosaic_df_loss$population
)

# Log transform for skew if needed
# mosaic_df_4326$log_housing_value <- log1p(mosaic_df_4326$housing_value)
# ------------------- Plot  -------------------------------

# --------- First, get bank average value in 2021 -----------------------------------------------
# Load the raster
bank_fill <- function(raster_path, polygon_sf, variable = "housing_units") {
  # Load raster and reproject to EPSG:4326
  raster_orig <- rast(raster_path)
  raster_4326 <- terra::project(raster_orig, "EPSG:4326")
  
  # Extract values with coordinates
  raster_df <- as.data.frame(raster_orig, xy = TRUE, na.rm = TRUE)
  colnames(raster_df)[3:5] <- c("housing_units", "housing_value", "population")
  
  # Reproject coordinates to EPSG:4326
  coords_4326 <- terra::project(cbind(raster_df$x, raster_df$y), crs(raster_orig), "EPSG:4326")
  
  # Combine into spatial points dataframe
  raster_df <- data.frame(
    x = coords_4326[,1],
    y = coords_4326[,2],
    housing_units = raster_df$housing_units,
    housing_value = raster_df$housing_value,
    population = raster_df$population
  )
  
  raster_points <- st_as_sf(raster_df, coords = c("x", "y"), crs = st_crs(polygon_sf))
  
  # Find points within polygon
  pts_in_poly <- raster_points[polygon_sf, , op = st_within]
  
  # Calculate mean of selected variable
  avg_val <- mean(pts_in_poly[[variable]], na.rm = TRUE)
  
  # Return the polygon with new variable added
  polygon_sf %>% mutate(!!variable := avg_val)
}

cypress_bank_filled <- bank_fill(
  raster_path = paste0("Service Area Mosaics/Big_Cypress_MB_Phase_I-V_mosaic.tif"),
  polygon_sf = cypress_bank_4326,
  variable = "housing_units"  # can also be "housing_value" or "population"
)

cypress_loss_filled <- bank_fill(
  raster_path = paste0("Service Area Mosaics/Loss Mosaics/Big_Cypress_MB_Phase_I-V_lossmosaic.tif"),
  polygon_sf = cypress_sa_4326,
  variable = "housing_units"  # can also be "housing_value" or "population"
)
cypress_loss_filled

# Load the raster
cypress_mosaic_sa <- rast(paste0("Service Area Mosaics/Big_Cypress_MB_Phase_I-V_mosaic.tif"))

# Extract values as a data frame (in native CRS, meters)
mosaic_df_sa <- as.data.frame(cypress_mosaic_sa, xy = TRUE, na.rm = TRUE)
colnames(mosaic_df_sa)[3:5] <- c("housing_units", "housing_value", "population")

# Reproject raster coordinates to EPSG:4326
coords_4326 <- terra::project(cbind(mosaic_df_sa$x, mosaic_df_sa$y), crs(cypress_mosaic_sa), "EPSG:4326")

# Create reprojected data frame for plotting
mosaic_df_sa <- data.frame(
  x = coords_4326[,1],
  y = coords_4326[,2],
  housing_units = mosaic_df_sa$housing_units,
  housing_value = mosaic_df_sa$housing_value,
  population = mosaic_df_sa$population
)


saplot_nbm = ggplot() +
  annotate("rect", xmin = xlim[1], xmax = xlim[2],
           ymin = ylim[1], ymax = ylim[2],
           fill = "lightblue3") +
  geom_sf(data = states, fill = "gray98", color = "gray70", linewidth = 0.2) +
  
  geom_sf(
    data  = cypress_sa_4326,
    aes(color = "Service Area"),
    fill = NA,
    linewidth = 1
  ) +
  
  geom_sf(
    data = cypress_bank_filled,
    aes(
      fill  = housing_units,
      color = "Bank, avg. value 2021"
    ),
    linewidth = 0.4
  ) +
  
  geom_tile(
    data   = mosaic_df_sa,
    aes(x = x, y = y, fill = housing_units),
    width  = lon_deg,
    height = lat_deg,
    alpha  = 0.8
  ) +
  
  geom_sf_label(
    data = places_fl %>% filter(NAME == "Cape Coral"),
    aes(label = NAME),
    stat = "sf_coordinates",
    nudge_y = 0.08,
    size = 3,
    label.size = 0.2,
    label.padding = unit(1, "mm")
  ) +
  
  geom_sf_label(
    data = places_fl %>% filter(NAME %in% c("Naples", "Fort Myers")),
    aes(label = NAME),
    stat = "sf_coordinates",
    size = 3,
    label.size = 0.2,
    label.padding = unit(1, "mm")
  ) +
  
  scale_fill_viridis_c(
    option = "viridis",
    name = "Downstream Housing Units",
    labels = scales::comma_format(),
    direction = 1
  ) +
  
  scale_color_manual(
    name = NULL,
    values = c(
      "Bank, avg. value 2021" = "red",
      "Service Area"          = "steelblue"
    ),
    guide = guide_legend(
      override.aes = list(
        fill = NA,
        linewidth = 1
      )
    )
  ) +
  
  guides(
    color = guide_legend(order = 1),
    fill  = guide_colorbar(order = 2)
  ) +
  
  coord_sf(xlim = xlim, ylim = ylim, expand = FALSE, crs = st_crs(target_crs)) +
  annotation_scale(location = "bl", width_hint = 0.25) +
  annotation_north_arrow(location = "tr", which_north = "true",
                         pad_x = unit(0.2, "in"), pad_y = unit(0.3, "in"),
                         style = north_arrow_fancy_orienteering) +
  theme_minimal(base_size = 13) +
  labs(
    title    = "Big Cypress Mitigation Bank Service Area",
    subtitle = "Remaining Wetlands (2021)",
    fill     = "Downstream Housing Units",
    x        = NULL,
    y        = NULL
  ) +
  theme(
    panel.grid       = element_blank(),
    legend.position  = "right",
    plot.title       = element_text(size = 13, hjust = 0.5),
    plot.subtitle    = element_text(size = 12, hjust = 0.5),
    legend.title     = element_text(size = 11),
    legend.text      = element_text(size = 11)
  )
saplot_nbm


# ---------------------- Same LOSS map WITHOUT basemap ----------------------

cell_res <- res(mosaic_raster_4326)  # returns a vector: c(xres, yres)
lon_deg <- cell_res[1]
lat_deg <- cell_res[2]

places_fl <- places(state = "FL", year = 2021)  # or latest year available

lossplot_nbm = ggplot() +
  # Draw ocean *first* with soft blue
  annotate(
    "rect",
    xmin = xlim[1], xmax = xlim[2],
    ymin = ylim[1], ymax = ylim[2],
    fill = "lightblue3"
  ) +
  
  # Then draw land on top
  geom_sf(data = states, fill = "gray98", color = "gray70", linewidth = 0.2) +
  
  # Service area outline (legend label "Service Area")
  geom_sf(
    data  = cypress_sa_4326,
    aes(color = "Service Area"),
    fill = NA,
    linewidth = 1
  ) +
  
  # Bank outline (legend label "Bank, avg. value 2021") + fill by housing_units
  geom_sf(
    data = cypress_bank_filled,
    aes(
      fill  = housing_units,          # <— Map the same fill you want to see in the legend
      color = "Bank, avg. value 2021" # <— The outline color goes to a separate legend item
    ),
    linewidth = 0.4
  )+
  
  # Downstream housing units (raster layer)
  geom_tile(
    data   = mosaic_df_4326_loss,  # Use the reprojected data
    aes(x = x, y = y, fill = housing_units),
    # Calculate appropriate width/height instead of using undefined variables
    width  = lat_deg,  # Example width in degrees (adjust as needed)
    height = lon_deg,  # Example height in degrees (adjust as needed)
    alpha  = 0.8
  )+
  
  geom_sf_label(
    data = places_fl %>% filter(NAME == "Cape Coral"),
    aes(label = NAME),
    stat = "sf_coordinates",
    nudge_y = 0.08,
    size = 3,
    label.size = 0.3,
    label.padding = unit(1, "mm")
  ) +
  
  # Naples and Fort Myers, no nudging
  geom_sf_label(
    data = places_fl %>% filter(NAME %in% c("Naples", "Fort Myers")),
    aes(label = NAME),
    stat = "sf_coordinates",
    size = 3,
    label.size = 0.2,
    label.padding = unit(1, "mm")
  )+
  
  # Fill scale for raster
  scale_fill_viridis_c(
    option = "viridis",
    name = "Downstream Housing Units",
    labels = scales::comma_format(),
    direction = 1
  ) +
  
  # Color scale for outlines
  scale_color_manual(
    name   = NULL,
    values = c(
      "Bank, avg. value 2021" = "red",
      "Service Area"          = "steelblue"
    ),
    # Override the default legend key appearance
    guide = guide_legend(
      override.aes = list(
        # Remove any fill for both keys in the color legend
        fill = NA,
        # Make both legend keys the same line width
        linewidth = 1
      )
    )
  )+
  guides(
    color = guide_legend(order = 1),
    fill  = guide_colorbar(order = 2)
  ) +
  
  coord_sf(xlim = xlim, ylim = ylim, expand = FALSE, crs = st_crs(target_crs)) +
  annotation_scale(location = "bl", width_hint = 0.25) +
  annotation_north_arrow(location = "tr", which_north = "true",
                         pad_x = unit(0.2, "in"), pad_y = unit(0.3, "in"),
                         style = north_arrow_fancy_orienteering) +
  
  theme_minimal(base_size = 13) +
  
  labs(
    title    = "Big Cypress Mitigation Bank Service Area",
    subtitle = "Lost Wetlands (1985–2021)",
    fill     = "Downstream Housing Units", 
    x        = NULL,
    y        = NULL
  ) +
  
  theme(
    panel.grid = element_blank(),  
    legend.position   = "right",
    plot.title        = element_text(size = 13, hjust = 0.5),
    plot.subtitle     = element_text(size = 12, hjust = 0.5),
    legend.title      = element_text(size = 11),
    legend.text       = element_text(size = 11)
  )

lossplot_nbm


# --------------------Before mapping, get bounding boxes -------------------------------

# 2) Get bounding box
bb <- st_bbox(cypress_sa_4326)
# Optionally pad a bit if you want
bb["xmin"] <- bb["xmin"] - 0.2
bb["xmax"] <- bb["xmax"] + 0.2
bb["ymin"] <- bb["ymin"] - 0.2
bb["ymax"] <- bb["ymax"] + 0.2

# 3) Use xlim and ylim
xlim <- c(bb["xmin"], bb["xmax"])
ylim <- c(bb["ymin"], bb["ymax"])


####################################################
###### FIGURE 5                               ######
####################################################

# Update file paths
ecdf_loss = readRDS(paste0(datadir, "ECDF Functions/SA ECDF Functions - Loss/ecdf_fn_Big_Cypress_MB_Phase_I-V.rds"))
ecdf_sa = readRDS(paste0(datadir, "ECDF Functions/SA ECDF Functions/ecdf_fn_Big_Cypress_MB_Phase_I-V.rds"))
ecdf_bank = readRDS(paste0(datadir, "ECDF Functions/Bank ECDF Functions/bank_ecdf_fn_Big_Cypress_MB_Phase_I-V.rds"))


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

# Plot: Only Loss vs Service Area if Bank is missing
df_hu = df_all %>% 
  filter(variable=="Housing Units")

library(scales)

fig4_hu = ggplot(df_hu, aes(x = value, y = ecdf, color = source)) +
  geom_line(linewidth = 1.2) +
  scale_x_continuous(
    limits = c(0, 120000), 
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
    x = "Downstream Housing Units",
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

fig4_hu

####################################################
###### FIGURE 6: DISTANCE DECAY FUNCTION PLOT ######
####################################################

distance_decay <- function(distance_km) {
  distance_km_values <- c(0, 10, 20, 30, 40)
  decay_values <- c(0.25, 1, 0.75, 1, 0)
  
  interp_values <- approx(
    x = distance_km_values,
    y = decay_values,
    xout = distance_km,
    rule = 2,
    yleft = 0,
    yright = 0,
    method = "linear"
  )$y
  return(interp_values)
}



distances <- seq(0, 50, length.out = 500)
decay_values <- distance_decay(distances)

# Plot the distance decay function
plot(
  distances, decay_values, type = "l", col = "blue", lwd = 2,
  xlab = "Distance to census tract centroid (km)", ylab = "Decay Value", main = "Distance Decay Function"
)

# Add defined points to the plot
points(
  x = c(0, 10, 20, 30, 40),
  y = c(0.25, 1, 0.75, 1, 0),
  col = "red", pch = 19)
