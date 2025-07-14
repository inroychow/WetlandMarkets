library(ggplot2)

##########################
######## FIGURE S1 #######
##########################

# ECDFs -----


results_long_normalized = read.csv(paste0(datadir,"\\normalized_ecdfs_loss.csv"))


results_hv_pop <- results_long_normalized %>%
  filter(comparison %in% c("normalized_loss_values"))

# Recalculate ecdf_summary for hv and pop 
ecdf_summary_hv_pop <- results_hv_pop %>%
  group_by(x, comparison) %>%
  summarize(mean_y = weighted.mean(y, w = n_obs), .groups = "drop")

a_hv_pop <- ggplot(
  results_hv_pop,
  aes(x = x, y = y, group = unit, lwd = n_obs)
) +
  geom_line(data = ecdf_summary_hv_pop, aes(x = x, y = mean_y, color = "Mean ECDF"), inherit.aes = FALSE, linewidth = 1) +
  
  
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

a_hv_pop


##########################
######## FIGURE S2 #######
##########################

# Violins

# FIND THE CODE FOR THIS

##########################
######## FIGURE S3 #######
##########################

library(tidyverse)
sa_dir <- "Extractions and Summaries/Extract Service Areas"
bank_dir <- "Extractions and Summaries/Extract Banks"
loss_dir <- "Wetland Loss/Extractions/Extract Service Area Loss"
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

ecdf_loss = readRDS("L:\\Wetland Flood Mitigation\\CONUS\\Wetland Loss\\Extractions\\Extract Service Area Loss\\SA Summaries Loss\\SA ECDF Functions - Loss\\ecdf_fn_Big_Cypress_MB_Phase_I-V.rds")
ecdf_sa = readRDS("L:\\Wetland Flood Mitigation\\CONUS\\Extractions\\Extract Service Areas 2021\\SA ECDF Functions\\ecdf_fn_Big_Cypress_MB_Phase_I-V.rds")
ecdf_bank = readRDS("L:\\Wetland Flood Mitigation\\CONUS\\Extractions\\Extract Banks 2021\\Bank ECDF Functions\\bank_ecdf_fn_Big_Cypress_MB_Phase_I-V.rds")


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

summary_loss = readRDS("CONUS/Summaries/loss_summary.rds")
summary_bank = readRDS("CONUS/Summaries/bank_summary.rds")
sa_summary = readRDS("CONUS/Summaries/sa_summary.rds")

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
                    ~ file.path("CONUS", "Extractions with Geom",
                                "Bank Extract with Geom",
                                glue("bank_extract_{.x}withgeom.rds")))
  
  sa_path   <- file.path("CONUS", "Extractions with Geom",
                         "SA Extract with Geom",
                         glue("sa_extract_{canon}.rds"))
  
  # loss lives in the same “sa_extract” folder according to your note
  loss_path <- file.path("CONUS", "Extractions with Geom",
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
                    ~ file.path("CONUS", "Extractions with Geom",
                                "Bank Extract with Geom",
                                glue("bank_extract_{.x}withgeom.rds")))
  
  sa_path   <- file.path("CONUS", "Extractions with Geom",
                         "SA Extract with Geom",
                         glue("sa_extract_{canon}.rds"))
  
  loss_path <- file.path("CONUS", "Extractions with Geom",
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

