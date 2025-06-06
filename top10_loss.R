
purrr::walk(c("future", "future.apply", "ggplot2", "tidyverse", "tibble", "progressr", 
              "data.table", "dplyr", "terra", "tigris", "sf", "scales", "ggspatial", 
              "rosm", "raster", "rnaturalearth", "tidycensus", "broom", "purrr", "glue"), library, character.only = TRUE)


setwd("L:\\Wetland Flood Mitigation")

sas= readRDS("Service Areas/ServiceAreas_agg.rds")
banks = readRDS("Footprints/footprints_and_buffers.rds")

# Raw data 

bank_path <- file.path("Extractions and Summaries", "Extract Banks",
                       glue("bank_extract_{bank}.rds"))
sa_path   <- file.path("Extractions and Summaries", "Extract Service Areas",
                       glue("sa_extract_{bank}.rds"))
loss_path <- file.path("Wetland Loss", "Extractions and Summaries",
                       glue("sa_extract_{bank}.rds"))

# Top 10 service areas with highest loss value
top10_banks = banks[banks$Name %in% c("Loxahatchee_MB", 
                                      "Little_Pine_Island_MB", 
                                      "Panther_Island_MB", 
                                      "Panther_Island_MB_-_Expansion", 
                                      "Florida_Wetlandsbank_at_Pembroke_Pines_MB", 
                                      "FP_AndL_Everglades_Phase_I_MB",
                                      "FP_AndL_Everglades_Phase_II_MB",
                                      "Big_Cypress_MB_Phase_I-V",
                                      "Big_Cypress_MB_Phase_VI",
                                      "Crosby_Island_Marsh_MB",
                                      "St_Johns_MB",
                                      "Katy_Prairie_Stream",
                                      "Crosby_Island_Marsh_MB",
                                      "Poa_Boy_MB"
),]

top10_sas <- sas[sas$ID %in%
                   c("Loxahatchee_MB", 
                     "Little_Pine_Island_MB", 
                     "Panther_Island_MB",
                     "Florida_Wetlandsbank_at_Pembroke_Pines_MB",
                     "FP_AndL_Everglades_Phase_II_MB",
                     "Big_Cypress_MB_Phase_I-V",
                     "Three_Lakes_Regional_MB_FDOT",
                     "Crosby_Island_Marsh_MB",
                     "St_Johns_MB",
                     "Katy_Prairie_Stream",
                     "Crosby_Island_Marsh_MB",
                     "Poa_Boy_MB"
                   ),]


# Get summary table

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

# Round large numbers for better readability
summary_sorted_fmt <- summary_sorted %>%
  mutate(across(contains("housing_units"), ~ scales::comma(.)))

# Output to LaTeX
kable(summary_sorted_fmt, format = "latex", booktabs = TRUE,
      caption = "Housing Units Summary for South Florida Wetland Areas",
      col.names = c("Bank / Area", "Sum", "Mean", "Max", "Type")) %>%
  kable_styling(latex_options = c("striped", "hold_position"))


############ T-test


# -------------- helper --------------
safe_read_rds_renamed <- function(path) {
  df <- tryCatch(readRDS(path), error = function(e) NULL)
  if (!is.null(df) && ncol(df) >= 4) {
    colnames(df)[2:4] <- c("housing_units", "housing_value", "population")
  }
  df
}


# -------------- main test function --------------
run_t_tests_for_bank <- function(bank_names) {
  canon <- bank_names[1]
  
  bank_paths <- map(bank_names,
                    ~ file.path("CONUS", "Extractions", "Extract Banks 2021",
                                glue("bank_extract_{.x}.rds")))
  
  bank_df <- bank_paths %>%
    map(safe_read_rds_renamed) %>%
    compact() %>%
    bind_rows()
  
  sa_path   <- file.path("CONUS", "Extractions", "Extract Service Areas 2021",
                         glue("sa_extract_{canon}.rds"))
  loss_path <- file.path("CONUS", "Wetland Loss", "Extractions",
                         "Extract Service Area Loss",
                         glue("sa_extract_{canon}.rds"))
  
  sa_df   <- safe_read_rds_renamed(sa_path)
  loss_df <- safe_read_rds_renamed(loss_path)
  
  if (any(vapply(list(bank_df, sa_df, loss_df), is.null, logical(1)))) {
    return(NULL)
  }
  
  metrics <- c("housing_units", "housing_value", "population")
  
  compare_group <- function(group_df, label) {
    map_dfr(metrics, \(metric) {
      if (all(is.na(group_df[[metric]]))) {
        return(NULL)
      }
      t <- t.test(loss_df[[metric]], group_df[[metric]])
      tibble(
        bank       = canon,
        metric     = metric,
        comparison = label,
        mean_loss  = mean(loss_df[[metric]],  na.rm = TRUE),
        mean_group = mean(group_df[[metric]], na.rm = TRUE),
        t_stat     = unname(t$statistic),
        df         = unname(t$parameter),
        p_value    = formatC(t$p.value, format = "e", digits = 10)
      )
    })
  }
  
  bind_rows(
    compare_group(bank_df, "bank"),
    compare_group(sa_df,   "service_area")
  )
}

# -------------- run everything --------------
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

all_ttests <- map_dfr(top10_names, run_t_tests_for_bank)

# Example: housing-unit rows only
hu_ttest <- all_ttests %>% filter(metric == "housing_units")
head(hu_ttest)

# hvpop_ttest = all_ttests %>% 
#   filter(metric %in% c("population", "housing_value"))
# head(hvpop_ttest)


#----------------------------------------

# ---------------------------------------

library(sf)
library(stringr)
library(purrr)

# Define directories
sa_dir <- "Extractions and Summaries/Extract Service Areas"
bank_dir <- "Extractions and Summaries/Extract Banks"
loss_dir <- "Wetland Loss/Extractions and Summaries/Extract Service Area Loss"


# Helper function to load matching files
load_rds_by_bankname <- function(dir, prefix, names_vec) {
  list.files(dir, pattern = paste0("^", prefix, ".*\\.rds$"), full.names = TRUE) |>
    keep(~ any(str_detect(.x, fixed(names_vec)))) |>
    map(readRDS)
}

# Load files
sa_rds_list   <- load_rds_by_bankname(sa_dir,   "sa_extract_",   top10_names)
bank_rds_list <- load_rds_by_bankname(bank_dir, "bank_extract_", top10_names)
loss_rds_list <- load_rds_by_bankname(loss_dir, "sa_extract_",   top10_names)

# Function to rename columns 2:4
rename_cols <- function(df) {
  names(df)[2:4] <- c("housing_units", "housing_value", "population")
  df
}
# Apply to each list
sa_rds_list   <- map(sa_rds_list, rename_cols)
bank_rds_list <- map(bank_combined, rename_cols)
loss_rds_list <- map(loss_rds_list, rename_cols)

bank_files <- list.files(bank_dir, pattern = "^bank_extract_.*\\.rds$", full.names = TRUE) |>
  keep(~ any(str_detect(.x, fixed(top10_names))))
bank_rds_list <- map(bank_files, readRDS)
bank_rds_list <- map(bank_combined, rename_cols)


sa_files <- list.files(sa_dir, pattern = "^sa_extract_.*\\.rds$", full.names = TRUE) |>
  keep(~ any(str_detect(.x, fixed(top10_names))))
sa_rds_list <- map(sa_files, readRDS)
sa_rds_list   <- map(sa_rds_list, rename_cols)

loss_files <- list.files(loss_dir, pattern = "^sa_extract_.*\\.rds$", full.names = TRUE) |>
  keep(~ any(str_detect(.x, fixed(top10_names))))
loss_rds_list <- map(loss_files, readRDS)
loss_rds_list <- map(loss_rds_list, rename_cols)


# Extract base name (remove path, prefix, suffix, and any "expansion/phase" details)
base_names_bank <- bank_files |>
  basename() |>
  str_remove("^bank_extract_") |>
  str_remove("\\.rds$") |>
  str_replace("_-.*", "") |>                    # remove expansion suffix
  str_replace("_Phase.*", "") |>                # remove phase suffix
  str_replace("_I-VI?$", "")                    # remove Roman numerals
base_names_bank

base_names_sa = sa_files |>
  basename() |>
  str_remove("^sa_extract_") |>
  str_remove("\\.rds$") |>
  str_replace("_-.*", "") |>                    # remove expansion suffix
  str_replace("_Phase.*", "") |>                # remove phase suffix
  str_replace("_I-VI?$", "")   


base_names_loss = loss_files |>
  basename() |>
  str_remove("^sa_extract_") |>
  str_remove("\\.rds$") |>
  str_replace("_-.*", "") |>                    # remove expansion suffix
  str_replace("_Phase.*", "") |>                # remove phase suffix
  str_replace("_I-VI?$", "") 

bank_combined <- split(bank_rds_list, base_names_bank) |>
  map(~ do.call(rbind, .x))




########


# Function to combine and tag data by source
combine_sources <- function(bank_df, sa_df, loss_df, bank_name) {
  bind_rows(
    mutate(bank_df, source = "bank"),
    mutate(sa_df,   source = "sa"),
    mutate(loss_df, source = "loss")
  ) |> mutate(bank_name = bank_name)
}

# Combine each set into one long dataframe
combined_all <- map2(seq_along(bank_rds_list), names(bank_rds_list), function(i, name) {
  combine_sources(
    bank_rds_list[[i]],
    sa_rds_list[[i]],
    loss_rds_list[[i]],
    bank_name = name
  )
}) |> bind_rows()

# Reshape for faceted plotting
long_df <- combined_all |>
  pivot_longer(cols = c(housing_units, housing_value, population),
               names_to = "variable", values_to = "value")

# Make ECDF plot (distribution curve) by bank, grouped by source
ggplot(long_df, aes(x = value, color = source)) +
  stat_ecdf() +
  facet_grid(variable ~ bank_name, scales = "free") +
  theme_minimal() +
  labs(title = "Distribution Curves by Bank and Source",
       x = "Value", y = "ECDF") +
  theme(axis.text.x = element_text(angle = 45, hjust = 1))

#######

# ------------------- DIST PLOT

library(tidyverse)

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
# Optional: sample large groups for plotting
# ----------------------------
max_rows <- 10000

plot_df_sampled <- plot_df %>%
  group_by(group) %>%
  mutate(n_rows = n()) %>%
  group_modify(~ if (.x$n_rows[1] > max_rows) slice_sample(.x, n = max_rows) else .x) %>%
  ungroup() %>%
  dplyr::select(-n_rows)

# ----------------------------
# Example plot
# ----------------------------
library(ggplot2)
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


###### SPAT AUTOCORR Dutilleul

library(sf)          
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
# run everything 
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

result <- map_dfr(FP, run_modified_t_tests_for_bank)
result
# housing-unit rows only
result_hu <- result %>% filter(metric == "housing_units")
head(result_hu)


result_hvpop = result %>% filter(metric %in% c("housing_value", "population"))


############# ONE TAILED



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
        loss_df  %>% select(x, y, val = !!sym(metric)) %>% mutate(grp = 1),
        group_df %>% select(x, y, val = !!sym(metric)) %>% mutate(grp = 0)
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
  c("FP_L_Everglades_Phase_I_MB", "FPL_Everglades_Phase_II_MB"))

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
