
###################################################################################################################
######                       Case Study: Big Cypress MB in Southwest Florida                       ################
######     has the most total downstream housing value of lost wetlands                            ################
###################################################################################################################

purrr::walk(c("future", "future.apply", "ggplot2", "tidyverse", "tibble", "progressr", 
              "data.table", "dplyr", "terra", "tigris", "sf", "scales", "ggspatial", 
              "rosm", "raster", "rnaturalearth"), library, character.only = TRUE)

setwd("L:\\Wetland Flood Mitigation")

loss_summary = readRDS("Wetland Loss/Extractions and Summaries/full_summary_loss.rds")
bank_summary = readRDS("Extractions and Summaries/Summaries Banks/full_summary_bank.rds")

sas= readRDS("Service Areas/ServiceAreas_agg.rds")
banks = readRDS("Bank Footprints/footprints_and_buffers.rds")


#Find which states intersect this SA ----------------------------------------
# lox = test %>% 
# filter(sa_id=="Loxahatchee_MB")

lox_banks = banks[banks$Name =="Loxahatchee_MB"]
lox_sa = sas[sas$ID == "Loxahatchee_MB",]

lox_mosaic_loss = rast("Wetland Loss/Service Area Mosaics Loss/Loxahatchee_MB_lossmosaic.tif")
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

# Create a new data frame with reprojected coordinates
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
  raster_path = "CONUS/Service Area Mosaics 2021/Loxahatchee_MB_mosaic.tif",
  polygon_sf = lox_bank_4326,
  variable = "housing_units"  # can also be "housing_value" or "population"
)

lox_loss_filled = bank_fill(
  raster_path = "CONUS/Wetland Loss/Service Area Mosaics Loss - 3.5/Loxahatchee_MB_lossmosaic.tif",
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
    width  = lon_deg,  # width in degrees 
    height = lat_deg,  # height in degrees 
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
lox_mosaic_sa <- rast("CONUS/Service Area Mosaics 2021/Loxahatchee_MB_mosaic.tif")

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

# ECDFs --------------------------
ecdf_loss = readRDS("L:\\Wetland Flood Mitigation\\CONUS\\Wetland Loss\\Extractions\\Extract Service Area Loss\\SA Summaries Loss\\SA ECDF Functions - Loss\\ecdf_fn_Loxahatchee_MB.rds")
ecdf_sa = readRDS("L:\\Wetland Flood Mitigation\\CONUS\\Extractions\\Extract Service Areas 2021\\SA ECDF Functions\\ecdf_fn_Loxahatchee_MB.rds")
ecdf_bank = readRDS("L:\\Wetland Flood Mitigation\\CONUS\\Extractions\\Extract Banks 2021\\Bank ECDF Functions\\bank_ecdf_fn_Loxahatchee_MB.rds")


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


# ---HV 

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
    title = "Loxahatchee Mitigation Bank",
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


# Pop 

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
    title = "Loxahatchee Mitigation Bank",
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


#########
# T-tests 

lox_sa_extract = readRDS("CONUS/Extractions and Summaries/Extract Service Areas/sa_extract_Loxahatchee_MB.rds")
lox_bank_extract = readRDS("CONUS/Extractions and Summaries/Extract Banks/bank_extract_Loxahatchee_MB.rds")
lox_loss_extract = readRDS("CONUS/Wetland Loss/Extractions and Summaries/sa_extract_Loxahatchee_MB.rds")

colnames(lox_sa_extract)[2:4] <- c("housing_units", "housing_value", "population")
colnames(lox_bank_extract)[2:4] <- c("housing_units", "housing_value", "population")
colnames(lox_loss_extract)[2:4] <- c("housing_units", "housing_value", "population")

lox_sa_extract$Source <- "Remaining Service Area (2021)"
lox_bank_extract$Source <- "Bank"
lox_loss_extract$Source <- "Service Area Loss (1985-2021)"

# Combine into one data frame
df <- bind_rows(lox_sa_extract, lox_bank_extract)
df <- bind_rows(df, lox_loss_extract)
library(dplyr)

lox_subset <- df %>%
  filter(Source %in% c("Bank", "Service Area Loss (1985-2021)")) %>%
  mutate(log_hu = log1p(housing_units))

t.test(housing_units ~ Source, data = lox_subset)
t.test(log_pop ~ Source, data = lox_subset)

