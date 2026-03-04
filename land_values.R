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

root_dir        <- "L:/Wetland Flood Mitigation"
sa_mosaic_dir   <- file.path(root_dir, "CONUS", "Service Area Mosaics 2021")
loss_mosaic_dir <- file.path(root_dir, "CONUS", "Wetland Loss", "Service Area Mosaics Loss")

nolte_landval <- rast(file.path(root_dir, "CONUS", "Land Values", "1 estimates", "places_fmv_vacant.tif"))

sas   <- readRDS(file.path(root_dir, "Service Areas", "ServiceAreas_agg.rds")) |> st_as_sf()
banks <- readRDS(file.path(root_dir, "Footprints", "footprints_and_buffers.rds")) |> st_as_sf()

bank_pool_map <- readRDS(file.path(root_dir, "ECDF Functions", "bank_pool_map.rds"))

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


