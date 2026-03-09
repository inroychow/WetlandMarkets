library(terra)
library(sf)
library(tigris)       # state polygons
library(progressr)
library(data.table)
library(future.apply)


 setwd("L:/Wetland Flood Mitigation")
states <- states(cb = TRUE)  
states <- states %>% 
  filter(!STUSPS%in%c("AK", "HI", "PR", "CNMI", "AS", "GU", "VI", "MP"))
plot(states$geometry)

sas = readRDS("Service Areas/ServiceAreas_agg.rds")
sas = st_as_sf(sas)

states = st_transform(states, st_crs(sas))


sas_states <- st_intersection(sas, states)

# Keep only relevant columns
sas_states <- sas_states %>%
  dplyr::select(ID = ID, state = STUSPS) # Replace `your_service_area_id_column`

sas_state_pairs <- sas_states %>%
  st_set_geometry(NULL) %>%
  distinct(ID, state)
# If a service area intersects multiple states, group them
sas_state_list <- sas_states %>%
  st_set_geometry(NULL) %>%
  distinct(ID, state)  # One row per ID-state pair
# Join back to original service area object if desired
sas <- left_join(sas, sas_state_list, by = c("ID" = "ID"))


library(dplyr)

# Step 1: Join `loss_summary` to `sas_state_list` to get the state(s) per service area

loss_summary = readRDS("CONUS/Summaries/loss_summary.rds")
loss_with_states <- loss_summary %>%
  left_join(sas_state_list, by = c("sa_id" = "ID"))


# Step 1: Count how many states each `sa_id` appears in
state_counts <- loss_with_states %>%
  distinct(sa_id, state) %>%
  count(sa_id)

# Step 2: Filter to `sa_id`s that appear in only one state
single_state_ids <- state_counts %>%
  filter(n == 1) %>%
  pull(sa_id)

# Step 3: Filter the original `loss_with_states` to keep only those service areas
loss_single_state <- loss_with_states %>%
  filter(sa_id %in% single_state_ids)


# Step 3: Summarize total loss value by state
onestate_summary <- loss_one_state %>%
  group_by(states) %>%
  summarise(total_sum_housing_units = sum(sum_housing_units, na.rm = TRUE),
            count_service_areas = n())



# MULTI STATES

setwd("L:/Wetland Flood Mitigation")


# ── paths ──────────────────────────────────────────────────────────────────
ras_dir   <- "CONUS/Wetland Loss/Service Area Mosaics Loss - 3.5"
out_dir   <- "CONUS/Wetland Loss/Extractions/Extract State"
if (!dir.exists(out_dir)) dir.create(out_dir, recursive = TRUE)


example = rast(paste0(ras_dir, "/AA_Shaw_lossmosaic.tif"))
crs = st_crs(example)
# ── read service-area polygons ─────────────────────────────────────────────
sas  <- readRDS("Service Areas/ServiceAreas_agg.rds") |> st_as_sf()

# ── state polygons (CONUS only) ────────────────────────────────────────────
states <- states(cb = TRUE) |>
  dplyr::filter(!STUSPS %in% c("AK","HI","PR","CNMI","AS","GU","VI","MP")) |>
  st_transform(st_crs(sas))

states = st_transform(states, crs)
# ── intersect once:  one row = one (SA ∩ State) piece ─────────────────────
sas_states <- st_intersection(
  sas  |> dplyr::select(sa_id = ID),          # keep only the SA id + geometry
  states |> dplyr::select(state = STUSPS)     # keep only the state abbrev
) |>                                          # result: one row per SA–state piece
  st_make_valid()                             # (robustness when clipping)
sas_states = st_transform(sas_states, crs)

# ── identify which SAs cross ≥2 states ─────────────────────────────────────
multi_state_ids <- sas_states |>
  st_set_geometry(NULL) |>
  dplyr::count(sa_id) |>
  dplyr::filter(n > 1) |>
  dplyr::pull(sa_id)

# ── list the TIFF mosaics you actually have on disk ────────────────────────
tifs   <- list.files(ras_dir, pattern = "_lossmosaic\\.tif$", full.names = FALSE)
tif_id <- sub("_lossmosaic\\.tif$", "", tifs)

# ── which multi-state SAs still need to be processed? ──────────────────────
done   <- sub("^sa_extract_(.*)\\.rds$", "\\1",
              list.files(out_dir, pattern = "^sa_extract_.*\\.rds$"))

to_do  <- intersect(multi_state_ids, tif_id) |>
  setdiff(done)

cat("✦ Multi-state SAs still to extract:", length(to_do), "\n")
if (length(to_do) == 0) {
  message("✅ Nothing left to do."); quit(save = "no")
}

# List of already processed files
done_files <- list.files(out_dir, pattern = "^sa_extract_.*\\.rds$")
done_ids <- sub("^sa_extract_(.*)\\.rds$", "\\1", done_files)

# All multi-state SA IDs
multi_state_ids <- unique(sas_states$sa_id[duplicated(sas_states$sa_id)])

states_to_do <- setdiff(states$STUSPS, sub("^state_extract_(.*)\\.rds$", "\\1",
                                           list.files(out_dir, pattern = "^state_extract_.*\\.rds$")))

future::plan(multisession, workers = 10)
handlers(global = TRUE)

with_progress({
  p <- progressor(along = states_to_do)
  
  future_lapply(states_to_do, function(st_abbrev) {
    p(st_abbrev)
    
    # All SA pieces in this state
    sa_pieces <- sas_states[sas_states$state == st_abbrev, ]
    
    # Group by SA
    sa_ids <- unique(sa_pieces$sa_id)
    
    result_list <- lapply(sa_ids, function(sa) {
      tif_path <- file.path(ras_dir, paste0(sa, "_lossmosaic.tif"))
      if (!file.exists(tif_path)) return(NULL)  # skip if raster doesn't exist
      
      r <- rast(tif_path)
      poly_part <- sa_pieces[sa_pieces$sa_id == sa, ]
      
      vals <- terra::extract(r, poly_part, fun = sum, na.rm = TRUE, ID = FALSE)
      if (is.null(colnames(vals)) || any(colnames(vals) == "")) {
        colnames(vals) <- paste0("sum_", names(r))
      }
      cbind(poly_part |> st_set_geometry(NULL), vals)
    })
    
    state_result <- rbindlist(result_list, use.names = TRUE, fill = TRUE)
    
    saveRDS(state_result, file = file.path(out_dir, paste0("state_extract_", st_abbrev, ".rds")))
    invisible(NULL)
  })
})

future::plan(sequential)
message("✅ State-level extraction complete.")

#### GET SUMMARY OF MULTI STATE

library(tidyverse)

folder_path <- "L:/Wetland Flood Mitigation/CONUS/Wetland Loss/Extractions/Extract State"
files <- list.files(folder_path, pattern = "^state_extract_.*", full.names = TRUE)

# Updated function to extract state code robustly
get_state <- function(filename) {
  base <- basename(filename)
  match <- str_match(base, "state_extract_([A-Z]{2})")[,2]
  return(match)
}

# Initialize summary list
summary_list <- list()

library(readr)      # fast readers
library(tibble)     # for tibble()
library(purrr)      # map() helpers
library(stringr)    # you’re already using str_match()

# ---- loop over files, one at a time ----
library(tidyverse)

summary_df <- map_dfr(files, function(fp) {
  st <- get_state(fp)
  
  # Try reading the RDS file safely
  dat <- tryCatch(readRDS(fp), error = function(e) NULL)
  
  # Skip if unreadable or not a data frame
  if (is.null(dat) || !is.data.frame(dat)) {
    message("Skipping file: ", basename(fp), " (unreadable or not a data frame)")
    return(tibble(state = st, housing_units_sum = NA_real_))
  }
  
  # Skip if column is missing
  if (!("housing_units" %in% colnames(dat))) {
    message("Skipping file: ", basename(fp), " (missing 'housing_units')")
    return(tibble(state = st, housing_units_sum = NA_real_))
  }
  
  # Remove rows with NA or zero housing_units
  valid_vals <- dat$housing_units[!is.na(dat$housing_units)]
  
  # Skip if nothing valid
  if (length(valid_vals) == 0) {
    message("Skipping file: ", basename(fp), " ('housing_units' all NA)")
    return(tibble(state = st, housing_units_sum = NA_real_))
  }
  
  # Sum and return
  hu_sum <- sum(valid_vals, na.rm = TRUE)
  tibble(state = st, housing_units_sum = hu_sum)
})

saveRDS(summary_df, "CONUS\\Wetland Loss\\Extractions\\Extract State\\multistate_summary.rds")


#### COMBINE 

multistate_summary = readRDS("CONUS\\Wetland Loss\\Extractions\\Extract State\\multistate_summary.rds")
multistate_summary<- multistate_summary %>%
  group_by(state) %>%
  summarise(
    housing_units_sum = sum(housing_units_sum, na.rm = TRUE),
    .groups = "drop"
  )

multistate_summary
onestate_summary


# Rename columns to prepare for join
multi <- multistate_summary %>%
  rename(total_multi = housing_units_sum)

one <- onestate_summary %>%
  rename(total_one = total_sum_housing_units) %>%
  rename(state = states)  # match column name to `multi`

# Full join and sum
combined_summary <- full_join(one, multi, by = "state") %>%
  mutate(
    total_all = rowSums(across(c(total_one, total_multi)), na.rm = TRUE)
  )

combined_summary

library(ggplot2)

# Optional: sort by total_all for better readability
combined_summary_sorted <- combined_summary %>%
  arrange(desc(total_all))

library(scales)

library(ggplot2)
library(scales)
library(dplyr)

# Sort descending
combined_summary_sorted <- combined_summary %>%
  arrange(total_all)

bar_plot <- ggplot(combined_summary_sorted, aes(x = reorder(state, total_all), y = total_all)) +
  geom_bar(stat = "identity", fill = "steelblue", width = 0.7) +
  labs(
    x = "State",
    y = "Total Housing Units (Billions)",
    title = "Total Housing Units by State\n(One-State + Multi-State Service Areas)"
  ) +
  coord_flip() +
  scale_y_continuous(
    labels = label_number(scale = 1e-9, suffix = "B", accuracy = 1),
    expand = expansion(mult = c(0, 0.05))  # adds space above bars
  ) +
  theme_bw(base_size = 14) +
  theme(
    axis.text.y = element_text(size = 12, margin = margin(r = 5)),
    axis.text.x = element_text(size = 12),
    axis.title = element_text(size = 16),
    plot.title = element_text(size = 16, hjust = 0.5),
    panel.grid.major.y = element_blank(),
    panel.grid.minor = element_blank()
  )

bar_plot


######## STATE UNION #############
# FULL COPY-PASTE SCRIPT – Option 1 (“union by state, no double counting”)
# -----------------------------------------------------------------------
# Adjust the three paths below to match your folder structure
library(sf)
library(terra)
library(tigris)
library(dplyr)
library(tibble)

# ── Set file paths ──
ras_dir  <- "CONUS/Wetland Loss/Service Area Mosaics Loss - 3.5"
sa_rds   <- "Service Areas/ServiceAreas_agg.rds"

# ── Load CONUS States and Service Areas ──
states <- states(cb = TRUE) |> 
  filter(!STUSPS %in% c("AK", "HI", "PR", "VI", "GU", "MP", "AS", "UM"))   # keep only CONUS

sas <- readRDS(sa_rds) |> st_as_sf()

# Reproject states to match service areas
states <- st_transform(states, st_crs(sas))

# Intersect service areas with each state
sas_states <- st_intersection(
  sas |> select(sa_id = ID),
  states |> select(state = STUSPS)
) |> st_make_valid()

# Build lookup table: state → list of SA IDs
sa_lookup <- sas_states |>
  st_set_geometry(NULL) |>
  group_by(state) |>
  summarise(sa_ids = list(unique(sa_id)), .groups = "drop")


# sas_per_state <- sas_states |>
#   st_drop_geometry() |>        # Drop geometry to speed up and simplify counting
#   distinct(sa_id, state) |>    # Make sure each SA–state pair is counted once
#   count(state, name = "n_sas") # Count how many unique SAs per state
# 
# print(sas_per_state)

# ── 2. helper: sum *new* cells only ---------------------------------------
sum_unique_cells <- function(r, processed_mask) {
  # r : SpatRaster for one SA, already cropped/masked to the state piece
  # processed_mask : SpatRaster (logical) indicating which cells were
  #                  already counted in previous SAs
  new_cells <- is.na(processed_mask) & !is.na(r)
  if (!any(new_cells[])) return(list(sum = 0, mask = processed_mask))  # nothing new
  
  # accumulate sum *only* from the new cells
  new_vals <- r
  new_vals[!new_cells] <- NA
  add_sum <- global(new_vals, fun = "sum", na.rm = TRUE)[1,1]
  
  # update the mask for next iteration
  processed_mask[new_cells] <- TRUE
  list(sum = add_sum, mask = processed_mask)
}

# ── 3. per-state loop (parallel) ------------------------------------------
plan(multisession, workers = 10)
handlers(global = TRUE)

state_totals <- with_progress({
  p <- progressor(along = nrow(sa_lookup))
  
  future_lapply(seq_len(nrow(sa_lookup)), function(i) {
    st_abbrev <- sa_lookup$state[i]
    sa_ids    <- sa_lookup$sa_ids[[i]]
    p(message = st_abbrev)
    
    # state polygon (single feature)
    state_poly  <- states[states$STUSPS == st_abbrev, ] |> st_make_valid()
    state_vect  <- vect(state_poly)                        # SpatVector
    
    # build an *empty mask* covering exactly that state’s extent
    # use the grid of the first raster we manage to open
    first_tif   <- NULL
    for (id in sa_ids) {
      fp <- file.path(ras_dir, paste0(id, "_lossmosaic.tif"))
      if (file.exists(fp)) { first_tif <- fp; break }
    }
    if (is.null(first_tif)) return(tibble(state = st_abbrev, total_loss = NA_real_))
    
    template_r   <- rast(first_tif)
    template_r   <- mask(template_r, state_vect)
    processed_ms <- classify(template_r, cbind(NA, NA))    # logical template mask
    
    # running total
    tot <- 0
    
    # loop over *this* state’s service areas
    for (sa in sa_ids) {
      fp <- file.path(ras_dir, paste0(sa, "_lossmosaic.tif"))
      if (!file.exists(fp)) next
      
      r   <- rast(fp)                     # open single SA raster
      r   <- mask(r, state_vect)          # keep only the part inside this state
      
      # skip if completely outside / all NA
      if (all(is.na(r[]))) { rm(r); next }
      
      # add the “new” cells only
      r <- terra::resample(r, processed_ms, method = "near")
      res <- sum_unique_cells(r, processed_ms)
      tot <- tot + res$sum
      processed_ms <- res$mask
      
      rm(r); gc()
    }
    
    tibble(state = st_abbrev, total_loss = tot)
  })
})

state_totals <- bind_rows(state_totals)
saveRDS(state_totals, out_rds)

plan(sequential)
message("✅ Finished – state_totals_nodup.csv written.")

top20 <- state_totals %>%
  slice_max(order_by = total_loss, n = 20)

bar_plot20 <- ggplot(top20, aes(x = reorder(state, total_loss), y = total_loss)) +
  geom_bar(stat = "identity", fill = "steelblue", width = 0.7) +
  labs(
    x = "State",
    y = "Total Housing Units",
    title = "Top 20 States by Total Housing Units\n(One-State + Multi-State Service Areas)"
  ) +
  coord_flip() +
  scale_y_continuous(
    expand = expansion(mult = c(0, 0.05))
  ) +
  theme_bw(base_size = 14) +
  theme(
    axis.text.y = element_text(size = 12, angle = 45, hjust = 1),
    axis.text.x = element_text(size = 12),
    axis.title = element_text(size = 16),
    plot.title = element_text(size = 16, hjust = 0.5),
    panel.grid.major.y = element_blank(),
    panel.grid.minor = element_blank()
  )

bar_plot20



###### MOSAICKING BY STATE

# ---------- 0.  libraries & parallel plan ----------
library(sf)           # vector I/O
library(terra)        # raster I/O
library(tigris)       # state polygons
library(future.apply)

options(tigris_use_cache = TRUE, tigris_class = "sf")
plan(multisession, workers = 5)   # adjust cores as you wish

# ---------- 1.  read layers ----------
sa_polygons  <- readRDS("Service Areas/ServiceAreas_agg.rds")  # contains field ID
states       <- states(cb = TRUE, year = 2021)                 # 50 states + DC
states       <- st_transform(states, st_crs(sa_polygons))      # match CRS
sa_polygons = st_as_sf(sa_polygons)
# ---------- 2.  build SA–state lookup ----------
# (SA polygons may cross state borders → intersection gives duplicates)
sa_state_lookup <- st_intersection(
  sa_polygons[, c("ID")],
  states[, c("STUSPS")]
) |>
  st_drop_geometry() |>
  split(~STUSPS)  # Makes it a list by state abbreviation

# make a vector of state abbreviations to iterate over
state_indices <- names(sa_state_lookup)

# ---------- 3.  helper for mosaicking to state ----------
output_state_dir <- "CONUS/Wetland Loss/State mosaics"
dir.create(output_state_dir, showWarnings = FALSE)

process_state <- function(state_abbr, lookup, sa_poly, out_dir) {
  library(terra)
  library(sf)
  message(Sys.time(), " | Starting ", state_abbr, " on PID: ", Sys.getpid())
  
  # path for final state mosaicr
  out_file <- file.path(out_dir, paste0(state_abbr, "_mosaic.tif"))
  if (file.exists(out_file)) {
    message("Skipping ", state_abbr, " (already processed).")
    return(out_file)
  }
  
  # ---- 3a. rasters to include ----
  sa_ids   <- lookup[[state_abbr]]$ID
  tif_paths <- file.path("CONUS/Wetland Loss/Service Area Mosaics Loss - 3.5", paste0(sa_ids, "_lossmosaic.tif"))
  tif_paths <- tif_paths[file.exists(tif_paths)]
  
  if (length(tif_paths) == 0) {
    warning("No SA mosaics found for ", state_abbr)
    return(NA)
  }
  
  message("Reading ", length(tif_paths), " SA rasters for ", state_abbr)
  ras_list <- lapply(tif_paths, rast)
  
  # ---- 3b. mosaic SAs together ----
  state_mosaic <- do.call(terra::mosaic, ras_list)
  
  # ---- 3c. crop/mask to state polygon ----
  state_poly <- vect(states[states$STUSPS == state_abbr, ])
  # terra object
  
  state_mosaic <- terra::crop(state_mosaic, state_poly)
  state_mosaic <- terra::mask(state_mosaic, state_poly)
  
  # ---- 3d. write result ----
  terra::writeRaster(state_mosaic, out_file, overwrite = FALSE)
  message("Finished ", state_abbr)
  out_file
}
head(state_indices)

plan(multisession, workers=2)
# ---------- 4.  run in parallel ----------
plan(multisession, workers=7)
results <- future_lapply(
  X   = state_indices,
  FUN = process_state,
  lookup = sa_state_lookup,
  sa_poly = vect(states),          # convert once for terra
  out_dir = output_state_dir,
  future.seed = TRUE
)

plan(sequential)


# # Test on a single state — e.g., Florida
# state_indices <- state_indices[!state_indices %in% c("FL", "AK", "WY", "CO")]
# tests =  state_indices[2:3]
# 
# results <- lapply(
#   X   = test_state,
#   FUN = process_state,
#   lookup = sa_state_lookup,        # ← already a list by state
#   sa_poly = vect(states),
#   out_dir = output_state_dir)



### EXTRACTION ###

library(terra)
library(progressr)
library(data.table)
library(future.apply)

setwd("L:/Wetland Flood Mitigation")

# Paths
raster_folder <- "CONUS/Wetland Loss/State mosaics"
output_folder <- "CONUS/Wetland Loss/Extractions/Extract Service Area Loss - State"

if (!dir.exists(output_folder)) {
  dir.create(output_folder, recursive = TRUE)
}


tif_files <- list.files(raster_folder, pattern = "_mosaic\\.tif$", full.names = FALSE)
tif_ids <- gsub("_mosaic\\.tif$", "", tif_files)

existing_files <- list.files(output_folder, pattern = "^loss_extract_.*\\.rds$", full.names = FALSE)
processed_state_ids <- gsub("^loss_extract_(.*)\\.rds$", "\\1", existing_files)

states_to_process <- state_indices[state_indices %in% tif_ids & !(state_indices %in% processed_state_ids) ]

cat("Total TIFF mosaics found:", length(tif_ids), "\n")
cat("Already processed:", length(processed_state_ids), "\n")
cat("Remaining to process:", nrow(states_to_process), "\n\n")

if (nrow(states_to_process) == 0) {
  cat("All mosaics with TIFFs have been processed. Exiting...\n")
  quit(save = "no")
}

### Parallel setup
future::plan(multisession, workers = 7)

progressr::with_progress({
  p <- progressr::progressor(along = states_to_process)
  
  results <- future_lapply(states_to_process, function(state_abbr) {
    library(terra)
    library(sf)
    
    p(sprintf("Processing state: %s", state_abbr))
    
    raster_path <- file.path(raster_folder, paste0(state_abbr, "_mosaic.tif"))
    output_path <- file.path(output_folder, paste0("loss_extract_", state_abbr, ".rds"))
    
    if (!file.exists(raster_path)) {
      warning("Raster not found for ", state_abbr)
      return(NULL)
    }
    
    # Load raster and polygon
    state_raster <- rast(raster_path)
    state_poly   <- states[states$STUSPS == state_abbr, ]
    state_vect   <- terra::vect(state_poly)
    
    # Extract values within state
    extracted_data <- terra::extract(
      x     = state_raster,
      y     = state_vect,
      cells = TRUE
    )
    
    extracted_data$state <- state_abbr
    
    # Save result
    saveRDS(extracted_data, output_path)
    
    return(NULL)
  })
})

message("✅ Extraction complete! Results saved to: ", output_folder)

future::plan(sequential)

### Summarize and Plot ###
output_folder <- "CONUS/Wetland Loss/Extractions/Extract Service Area Loss - State"        # adjust if needed
file_list <- list.files(output_folder,
                        pattern = "^loss_extract_.*\\.rds$",
                        full.names = TRUE)

library(data.table)

dir.create("CONUS/Wetland Loss/Extractions/State summaries", showWarnings = FALSE)
plan(multisession, workers = 7)


# ---------------------------------------------------------------------------
# 1.  Process each RDS in parallel
# ---------------------------------------------------------------------------
state_summaries <- future_lapply(file_list, function(f) {
  
  dt <- readRDS(f)
  setDT(dt)
  
  # Rename the three duplicated 'lyr1' columns
  colnames(dt)[2:4] <- c("housing_units", "housing_value", "population")
  
  # Summarise by state
  ss <- dt[, .(
    cell_count        = .N,
    total_housing_units = sum(housing_units, na.rm = TRUE),
    mean_housing_units  = mean(housing_units, na.rm = TRUE),
    total_housing_value = sum(housing_value, na.rm = TRUE),
    mean_housing_value  = mean(housing_value, na.rm = TRUE),
    total_population    = sum(population, na.rm = TRUE),
    mean_population     = mean(population, na.rm = TRUE)
  ), by = state]
  
  # Save per-file summary
  out_name <- sub("\\.rds$", "_summaryupdate.rds", basename(f))
  saveRDS(ss, file.path(output_folder, "State summaries", out_name))
  
  ss  # Return to rbind later
})

# Georgia extract won't load in .rds, so get summary stats straight from mosaic:

library(terra)
library(data.table)

# Load Georgia raster
ga_raster <- rast("CONUS/Wetland Loss/Mosaics/housing_units_GA.tif")  # adjust path if needed

# Extract each layer by index
hu_layer   <- ga_raster[[1]]  # housing units
hv_layer   <- ga_raster[[2]]  # housing value
pop_layer  <- ga_raster[[3]]  # population

# Compute statistics
cell_count          <- global(!is.na(hu_layer), fun = "sum")$sum
total_housing_units <- global(hu_layer, fun = "sum", na.rm = TRUE)$sum
mean_housing_units  <- global(hu_layer, fun = "mean", na.rm = TRUE)$mean

total_housing_value <- global(hv_layer, fun = "sum", na.rm = TRUE)$sum
mean_housing_value  <- global(hv_layer, fun = "mean", na.rm = TRUE)$mean

total_population    <- global(pop_layer, fun = "sum", na.rm = TRUE)$sum
mean_population     <- global(pop_layer, fun = "mean", na.rm = TRUE)$mean

# Assemble into a single-row data.table
ga_summary <- data.table(
  state                = "GA",
  cell_count           = cell_count,
  total_housing_units  = total_housing_units,
  mean_housing_units   = mean_housing_units,
  total_housing_value  = total_housing_value,
  mean_housing_value   = mean_housing_value,
  total_population     = total_population,
  mean_population      = mean_population
)

print(ga_summary)
saveRDS(ga_summary, "CONUS\\Wetland Loss\\Extractions\\Extract Service Area Loss - State\\State summaries\\loss_extract_GA_summaryupdate.rds")
# ---------------------------


summary_dir <- file.path(output_folder, "State summaries")
# List all summary files
summary_files <- list.files(summary_dir, pattern = "_summaryupdate\\.rds$", full.names = TRUE)
all_summaries <- rbindlist(lapply(summary_files, readRDS))

# Plot 
library(scales)

all_summaries %>%
  arrange(desc(total_housing_units)) %>%
  slice_head(n = 20) %>%
  ggplot(aes(x = reorder(state, total_housing_units), y = total_housing_units)) +
  geom_col(fill = "steelblue") +
  labs(
    x = "State",
    y = "Downstream Housing Units"
  ) +
  coord_flip() +
  scale_y_continuous(
    labels = label_number(scale = 1e-9, suffix = "T")
  ) +
  theme_minimal(base_size = 14) +
  theme(
    panel.grid.major.y = element_blank(),
    plot.title = element_text(face = "bold")
  )

library(dplyr)
library(ggplot2)
library(scales)

all_summaries %>%
  arrange(desc(total_housing_units)) %>%
  slice_head(n = 20) %>%
  ggplot(aes(x = reorder(state, total_housing_units), y = total_housing_units)) +
  geom_col(fill = "steelblue") +
  labs(
    title = "Top 20 States by Total Housing Units (Log Scale)",
    x = "State",
    y = "Housing Units (log scale)"
  ) +
  coord_flip() +
  scale_y_log10(
    labels = label_number(scale = 1e-9)
  ) +
  theme_minimal(base_size = 14) +
  theme(
    panel.grid.major.y = element_blank(),
    plot.title = element_text(face = "bold")
  )


library(data.table)

setDT(all_summaries)

# Filter Florida
fl <- all_summaries[state == "FL"]
tx <- all_summaries[state == "TX"]

# Total across all states
total_units <- sum(all_summaries$total_housing_units, na.rm = TRUE)
total_cells = sum(all_summaries$cell_count, na.rm = TRUE)
# Percentages
YY <- (tx$total_housing_units / total_units) * 100
XX = (fl$cell_count / total_cells) *100

all_summaries[order(-cell_count)][1:5, .(state, cell_count)]


# Round and print
cat(sprintf("Florida accounts for %.1f%% of total lost flood-protective value (TX)\n", YY))
ca_raster <- rast("CONUS/Wetland Loss/State mosaics/CA_mosaic.tif")  # adjust path if needed
plot(ca_raster[[1]])

# -----------------------

# CHANGE IN HOUSING UNITS #


#########

library(tidycensus)
library(tigris)
library(dplyr)
library(purrr)
library(sf)
library(tmap)

options(tigris_use_cache = TRUE)

conus_states <- fips_codes %>%
  filter(!state %in% c("AK", "HI", "PR", "AS", "GU", "MP", "VI", "UM")) %>%
  distinct(state) %>%
  pull(state)
# Wrap safely
safe_get_decennial <- safely(function(state_abbr) {
  get_decennial(
    geography = "tract",
    variables = "H001001",
    state = state_abbr,
    year = 2000,
    geometry = FALSE
  ) %>%
    rename(housing_units_2000 = value)
})

hu_2000_list <- map(conus_states, ~ safe_get_decennial(.x)$result)
# Remove NULL results
hu_2000_conus <- bind_rows(hu_2000_list[!map_lgl(hu_2000_list, is.null)])



# Wrap safely
safe_get_acs <- safely(function(state_abbr) {
  get_acs(
    geography = "tract",
    variables = "B25001_001",
    state = state_abbr,
    year = 2021,
    geometry = FALSE
  ) %>%
    rename(housing_units_2021 = estimate)
})

hu_2021_list <- map(conus_states, ~ safe_get_acs(.x)$result)
# Remove NULL results
hu_2021_conus <- bind_rows(hu_2021_list[!map_lgl(hu_2021_list, is.null)])


state_2000 <- hu_2000_conus %>%
  mutate(state = substr(GEOID, 1, 2)) %>%
  left_join(state_lookup, by = "state") %>%
  group_by(state_name) %>%
  summarize(total_2000 = sum(housing_units_2000, na.rm = TRUE))

state_2021 <- hu_2021_conus %>%
  mutate(state = substr(GEOID, 1, 2)) %>%
  left_join(state_lookup, by = "state") %>%
  group_by(state_name) %>%
  summarize(total_2021 = sum(housing_units_2021, na.rm = TRUE))

state_summary_full <- left_join(state_2000, state_2021, by = "state_name") %>%
  mutate(
    total_change = total_2021 - total_2000,
    pct_change = 100 * total_change / total_2000
  ) %>%
  arrange(desc(pct_change))
state_summary_full
