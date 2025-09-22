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

######################################
######          FIGURE 2        ######
######                          ######
######################################

suppressPackageStartupMessages({
  library(tidyverse)
  library(ggplot2)
  library(scales)
  library(patchwork)
  library(vctrs)
})

compare_remaining_vs_lost <- function(datadir, sample_size = 10000, seed = 123, bank_filter = NULL) {
  suppressPackageStartupMessages({
    library(tidyverse)
    library(scales)
    library(vctrs)
  })
  
  message("Scanning directory for bank ECDF files …")
  bank_files <- list.files(file.path(datadir, "SA ECDF Functions"), full.names = FALSE)
  bank_names <- sub("^ecdf_fn_(.*)\\.rds$", "\\1", bank_files)
  
  
  # ✅ Filter bank_names if filter is provided
  if (!is.null(bank_filter)) {
    bank_names <- bank_names[bank_names %in% bank_filter]
  }
  
  total_banks <- length(bank_names)
  if (total_banks == 0) stop("No ECDF files found for specified banks in directory: ", datadir)
  
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

# -------------

###########################################################
############# FIGURE 4 ECDFs Housing Units ################
###########################################################
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
bank_names <- sub("^ecdf_fn_(.*)\\.rds$", "\\1", name_files)

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

################# 
### Figure 4 ####
#################
# ECDF Lox and Cypress HU



# ECDFs --------------------------
ecdf_loss = readRDS(paste0(datadir,"Extractions and Summaries/Loss ECDF Functions/ecdf_fn_Loxahatchee_MB.rds"))
ecdf_sa = readRDS(paste0(datadir, "Extractions and Summaries/SA ECDF Functions\\ecdf_fn_Loxahatchee_MB.rds"))
ecdf_bank = readRDS(paste0(datadir, "Extractions and Summaries/Bank ECDF Functions\\bank_ecdf_fn_Loxahatchee_MB.rds"))


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

ecdf_loss = readRDS("Extractions and Summaries/Loss ECDF Functions/ecdf_fn_Big_Cypress_MB_Phase_I-V.rds")
ecdf_sa = readRDS("Extractions and Summaries/SA ECDF Functions/ecdf_fn_Big_Cypress_MB_Phase_I-V.rds")
ecdf_bank = readRDS("Extractions and Summaries/Bank ECDF Functions/bank_ecdf_fn_Big_Cypress_MB_Phase_I-V.rds")


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

loss_summary = readRDS("CONUS/Wetland Loss/Extractions/Extract Service Area Loss/SA Summaries Loss/sa_summary_loss.rds")
bank_summary = readRDS("Extractions/Extract Banks 2021/bank_summary_2021.rds")

sas= readRDS("Service Areas/ServiceAreas_agg.rds")
banks = readRDS("Footprints/footprints_and_buffers.rds")


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

test = readRDS(paste0(datadir,"Extractions and Summaries/sa_summary_loss.rds"))
test_bank = readRDS(paste0(datadir,"Extractions and Summaries/bank_summary.rds"))

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
ecdf_loss = readRDS(paste0(datadir, "Extractions and Summaries/Loss ECDF Functions/ecdf_fn_Big_Cypress_MB_Phase_I-V.rds"))
ecdf_sa = readRDS(paste0(datadir, "Extractions and Summaries/SA ECDF Functions/ecdf_fn_Big_Cypress_MB_Phase_I-V.rds"))
ecdf_bank = readRDS(paste0(datadir, "Extractions and Summaries/Bank ECDF Functions/bank_ecdf_fn_Big_Cypress_MB_Phase_I-V.rds"))


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
