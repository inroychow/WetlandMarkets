#####################################
setwd("C:\\Users\\indumati\\Box\\Paper2_final")

# Summary bank
combined_df = readRDS("Extractions and Summaries\\Extract Banks\\all_bank_extractions.rds")
library(dplyr)

# Load your data
combined_df <- readRDS("Extractions and Summaries/Extract Banks/all_bank_extractions.rds")
poolmap <- readRDS("Bank Footprints/pool_map.rds")

# Add pool information to extracted data
combined_df_pooled <- combined_df %>%
  left_join(poolmap %>% select(bank_name, pool_id), 
            by = "bank_name") %>%
  mutate(pool_id = coalesce(pool_id, paste0("Solo_", bank_name)))

# Aggregate by pool (treating all banks in same pool as one unit)
bank_pooled_summary <- combined_df_pooled %>%
  group_by(pool_id) %>%
  summarise(
    banks_in_pool = paste(unique(bank_name), collapse = "; "),
    n_banks = n_distinct(bank_name),
    sum_housing_units = sum(housing_units, na.rm = TRUE),
    mean_housing_units = mean(housing_units, na.rm = TRUE),
    median_housing_units = median(housing_units, na.rm = TRUE),
    min_housing_units = min(housing_units, na.rm = TRUE),
    max_housing_units = max(housing_units, na.rm = TRUE),
    sd_housing_units = sd(housing_units, na.rm = TRUE),
    q25_housing_units = quantile(housing_units, 0.25, na.rm = TRUE),
    q75_housing_units = quantile(housing_units, 0.75, na.rm = TRUE),
    
    sum_housing_value = sum(housing_value, na.rm = TRUE),
    mean_housing_value = mean(housing_value, na.rm = TRUE),
    median_housing_value = median(housing_value, na.rm = TRUE),
    min_housing_value = min(housing_value, na.rm = TRUE),
    max_housing_value = max(housing_value, na.rm = TRUE),
    sd_housing_value = sd(housing_value, na.rm = TRUE),
    q25_housing_value = quantile(housing_value, 0.25, na.rm = TRUE),
    q75_housing_value = quantile(housing_value, 0.75, na.rm = TRUE),
    
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
  ) %>%
  mutate(pool_id = gsub("^Solo_", "", pool_id))  # Remove Solo_ prefix

# Save pooled summary
# saveRDS(bank_pooled_summary, "Extractions and Summaries/bank_summary_pool.rds")


# Loss summary pooled -----------------
oldloss = readRDS("Extractions and Summaries/loss_summary.rds")

# Add pool information to loss_summary
loss_summary_pooled <- oldloss %>%
  left_join(poolmap %>% select(bank_name, pool_id), 
            by = c("sa_id" = "bank_name")) %>%
  mutate(pool_id = coalesce(pool_id, paste0("Solo_", sa_id)),
         pool_id = gsub("^Solo_", "", pool_id))

# Select one SA per pool
loss_summary_pooled <- loss_summary_pooled %>%
  group_by(pool_id) %>%
  slice(1) %>%  # Takes the first row per pool
  ungroup()


# Load both datasets
bank_pooled_summary <- readRDS("Extractions and Summaries/bank_summary_pool.rds")

# Get the pool_ids that exist in pooled_summary
pools_in_bank_summary <- unique(pooled_summary$pool_id)

# Filter loss_summary_pooled to only those pools
loss_summary_matched <- loss_summary_pooled %>%
  filter(pool_id %in% pools_in_bank_summary)

# saveRDS(loss_summary_matched, "Extractions and Summaries/loss_summary_pool.rds")
