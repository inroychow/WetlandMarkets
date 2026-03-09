# Cypress has the most lost wetland cells and the highest value lost in SA compared to bank.------------------------------------

purrr::walk(c("future", "future.apply", "ggplot2", "tidyverse", "tibble", "progressr", 
              "data.table", "dplyr", "terra", "tigris", "sf", "scales", "ggspatial", 
              "rosm", "raster", "rnaturalearth"), library, character.only = TRUE)

setwd("L:\\Wetland Flood Mitigation")

test = readRDS("CONUS/Wetland Loss/Extractions/Extract Service Area Loss/SA Summaries Loss/sa_summary_loss.rds")
test_bank = readRDS("Extractions/Extract Banks 2021/bank_summary_2021.rds")

sas= readRDS("Service Areas/ServiceAreas_agg.rds")
banks = readRDS("Footprints/footprints_and_buffers.rds")


#Find which states intersect this SA ----------------------------------------
# cypress = test %>% 
  # filter(sa_id=="Big_Cypress_MB_Phase_I-V")

cypress_banks = banks[banks$Name %in% c("Big_Cypress_MB_Phase_I-V", "Big_Cypress_MB_Phase_VI")]
cypress_banks = terra::union(cypress_banks)
cypress_sa = sas[sas$ID == "Big_Cypress_MB_Phase_I-V",]

cypress_mosaic_loss = rast("CONUS/Wetland Loss/Service Area Mosaics Loss - 3.5/Big_Cypress_MB_Phase_I-V_lossmosaic.tif")

states <- ne_states(country = "United States of America", returnclass = "sf")
states_albers <- st_transform(states, crs(cypress_sa))
states_vect <- vect(states_albers)
intersected <- terra::intersect(states_vect, cypress_sa)
intersected$name

fl_sf <- subset(states_albers, name %in% "Florida")

mosaic_df_loss <- as.data.frame(cypress_mosaic_loss, xy = TRUE, na.rm = TRUE)
colnames(mosaic_df_loss[3:5]) <- c("housing_units", "housing_value", "population")

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

# ------------------- Plot  -------------------------------

# --------- First, get bank average value in 2021 -----------------------------------------------
# Load the raster
bank_fill <- function(raster_path, polygon_sf, variable = "housing_value") {
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
  raster_path = "CONUS/Service Area Mosaics 2021/Big_Cypress_MB_Phase_I-V_mosaic.tif",
  polygon_sf = cypress_bank_4326,
  variable = "housing_value"  # can also be "housing_value" or "population"
)


# Lost wetlands, downstream Housing units -------------------------------------------------------------------------------

lossplot_hu = ggplot() +
  #geom_sf(data = states, fill = NA, color = "gray40", linewidth = 0.5) +
  annotation_map_tile(type = "cartolight", zoomin = -1) +
  
  # Service area outline (legend label "Service Area")
  geom_sf(
    data  = cypress_sf_4326,
    aes(color = "Service Area"),
    fill = NA,
    linewidth = 1
  ) +
  
  # Bank outline (legend label "Bank, avg. value 2021") + fill by housing_units
  geom_sf(
    data = cypress_bank_filled,
    aes(
      fill  = housing_value,          # <— Map the same fill you want to see in the legend
      color = "Bank, avg. value 2021" # <— The outline color goes to a separate legend item
    ),
    linewidth = 0.4
  )+

  # Downstream housing units (raster layer)
  geom_tile(
    data   = mosaic_df_loss,
    aes(x = x, y = y, fill = housing_value),
    width  = lon_deg,
    height = lat_deg,
    alpha  = 0.8
  ) +
  geom_sf_label(data = cities, aes(label = NAME), size = 3, label.size = 0.2, label.padding = unit(1, "mm"))+
  
  # Fill scale for raster
  scale_fill_viridis_c(
    option = "viridis",
    name = "Downstream Total Housing Value",
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
        # Or you can set shape, size, etc. if you want them as points or bigger squares, etc.
        
      )
    )
  )+
  
  coord_sf(xlim = xlim, ylim = ylim, expand = FALSE, crs = st_crs(target_crs)) +
  annotation_scale(location = "bl", width_hint = 0.25) +
  annotation_north_arrow(location = "tr", which_north = "true",
                         pad_x = unit(0.2, "in"), pad_y = unit(0.3, "in"),
                         style = north_arrow_fancy_orienteering) +
  
  theme_minimal(base_size = 13) +
  
  labs(
    title    = "Big Cypress Mitigation Bank Service Area",
    subtitle = "Lost Wetlands (1985–2021)",
    fill     = "Downstream Total Housing Value", 
    x        = NULL,
    y        = NULL
  ) +
  
  theme(
    legend.position   = "right",
    plot.title        = element_text(size = 13, hjust = 0.5),
    plot.subtitle     = element_text(size = 12, hjust = 0.5),
    legend.title      = element_text(size = 11),
    legend.text       = element_text(size = 11)
  )

lossplot_hu


# ---------------------- Same LOSS map WITHOUT basemap ----------------------
places_fl <- places(state = "FL", year = 2021)  # or latest year available

lossplot_nbm = ggplot() +
  # Draw ocean *first* with soft blue
  annotate(
    "rect",
    xmin = xlim[1], xmax = xlim[2],
    ymin = ylim[1], ymax = ylim[2],
    fill = "aliceblue"
  ) +
  
  # Then draw land on top
  geom_sf(data = states, fill = "gray98", color = "gray70", linewidth = 0.2) +
  
  # Service area outline (legend label "Service Area")
  geom_sf(
    data  = cypress_sf_4326,
    aes(color = "Service Area"),
    fill = NA,
    linewidth = 1
  ) +
  
  # Bank outline (legend label "Bank, avg. value 2021") + fill by housing_units
  geom_sf(
    data = cypress_bank_filled,
    aes(
      fill  = population,          # <— Map the same fill you want to see in the legend
      color = "Bank, avg. value 2021" # <— The outline color goes to a separate legend item
    ),
    linewidth = 0.4
  )+
  
  # Downstream housing units (raster layer)
  geom_tile(
    data   = mosaic_df_loss,
    aes(x = x, y = y, fill = population),
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
    name = "Downstream Population",
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
  
  coord_sf(xlim = xlim, ylim = ylim, expand = FALSE, crs = st_crs(target_crs)) +
  annotation_scale(location = "bl", width_hint = 0.25) +
  annotation_north_arrow(location = "tr", which_north = "true",
                         pad_x = unit(0.2, "in"), pad_y = unit(0.3, "in"),
                         style = north_arrow_fancy_orienteering) +
  
  theme_minimal(base_size = 13) +
  
  labs(
    title    = "Big Cypress Mitigation Bank Service Area",
    subtitle = "Lost Wetlands (1985–2021)",
    fill     = "Downstream Population", 
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
cypress_mosaic_sa <- rast("CONUS/Service Area Mosaics 2021/Big_Cypress_MB_Phase_I-V_mosaic.tif")

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


res_vals <- res(cypress_mosaic_4326)  # returns degrees
lon_deg <- res_vals[1]
lat_deg <- res_vals[2]

#  set extent from a relevant sf layer
zoom_ext_4326 <- st_bbox(cypress_sf_4326)
bb <- zoom_ext_4326
bb["xmin"] <- bb["xmin"] - 0.2
bb["xmax"] <- bb["xmax"] + 0.2
bb["ymin"] <- bb["ymin"] - 0.2
bb["ymax"] <- bb["ymax"] + 0.2
xlim <- c(bb["xmin"], bb["xmax"])
ylim <- c(bb["ymin"], bb["ymax"])

# Plot
saplot_hu = ggplot() +
  # Basemap
  annotation_map_tile(type = "cartolight", zoomin = -1) +
  
  geom_tile(
    data = mosaic_df_sa,
    aes(x = x, y = y, fill = housing_units),
    width  = lon_deg * 1.001,
    height = lat_deg * 1.001,
    color  = NA,
    size   = 0,
    alpha  = 1.0
  ) +
  
  # Service area outline
  geom_sf(
    data  = cypress_sf_4326,
    aes(color = "Service Area"),
    fill      = NA,
    linewidth = 1
  ) +
  
  # Bank outline + fill by housing_units
  geom_sf(
    data = cypress_bank_filled,
    aes(
      fill  = housing_units,
      color = "Bank, avg. value 2021"
    ),
    linewidth = 0.4
  ) +
  
  # Viridis fill scale
  scale_fill_viridis_c(
    option   = "viridis",
    name     = "Downstream Housing Units",
    labels   = scales::comma_format(),
    direction = 1
  ) +
  
  # Manual color scale for outlines
  scale_color_manual(
    name   = NULL,
    values = c(
      "Bank, avg. value 2021" = "red",
      "Service Area"          = "steelblue"
    ),
    guide = guide_legend(
      override.aes = list(
        fill      = NA,
        linewidth = 1
      ),
      order = 2  # This sets the order - higher number means lower position
    )
  ) +
  
  # Coordinate limits
  coord_sf(
    xlim   = xlim,
    ylim   = ylim,
    expand = FALSE,
    crs    = st_crs(4326)
  ) +
  
  # Scale bar and north arrow
  annotation_scale(location = "bl", width_hint = 0.25) +
  annotation_north_arrow(
    location = "tr",
    which_north = "true",
    pad_x = unit(0.2, "in"),
    pad_y = unit(0.3, "in"),
    style = north_arrow_fancy_orienteering
  ) +
  
  # Theme + labels
  theme_minimal(base_size = 13) +
  labs(
    title    = "Remaining Wetlands (2021) in Big Cypress Mitigation Bank Service Area",
    subtitle = "Flood protection value of remaining 2021 wetlands, based on housing units downstream",
    fill     = "Downstream Housing Units",
    x        = NULL,
    y        = NULL
  ) +
  theme(
    legend.position = "right",
    plot.title      = element_text(size = 13, hjust = 0.5),
    plot.subtitle   = element_text(size = 12, hjust = 0.5),
    legend.title    = element_text(size = 11),
    legend.text     = element_text(size = 11)
  ) +
  guides(
    fill = guide_colorbar(order = 1)  # This sets the housing units legend order - lower number means higher position
  )

saplot_hu

# Reproject the loss raster
cypress_mosaic_4326_loss <- terra::project(cypress_mosaic_loss, "EPSG:4326")

# Extract values + coords
mosaic_df_loss <- as.data.frame(cypress_mosaic_4326_loss, xy = TRUE, na.rm = TRUE)
colnames(mosaic_df_loss)[3:5] <- c("housing_units", "housing_value", "population")

# Raster resolution in degrees (use the reprojected loss raster!)
res_vals_loss <- res(cypress_mosaic_4326_loss)
lon_deg <- res_vals_loss[1]
lat_deg <- res_vals_loss[2]

lossplot_hu = ggplot() +
  #geom_sf(data = states, fill = NA, color = "gray40", linewidth = 0.5) +
  annotation_map_tile(type = "cartolight", zoomin = -1) +
  
  # Service area outline (legend label "Service Area")
  geom_sf(
    data  = cypress_sf_4326,
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
    data   = mosaic_df_loss,
    aes(x = x, y = y, fill = housing_units),
    width  = lon_deg,
    height = lat_deg,
    alpha  = 0.8
  ) +
  geom_sf_label(data = cities, aes(label = NAME), size = 3, label.size = 0.2, label.padding = unit(1, "mm"))+
  
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
        # Or you can set shape, size, etc. if you want them as points or bigger squares, etc.
        
      )
    )
  )+
  
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
    legend.position   = "right",
    plot.title        = element_text(size = 13, hjust = 0.5),
    plot.subtitle     = element_text(size = 12, hjust = 0.5),
    legend.title      = element_text(size = 11),
    legend.text       = element_text(size = 11)
  )

lossplot_hu
#  --------------------------------Plot remaining SAs mosaic -----------------------------------------------------

# Load the raster
cypress_mosaic_sa <- rast("CONUS/Service Area Mosaics 2021/Big_Cypress_MB_Phase_I-V_mosaic.tif")

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


res_vals <- res(cypress_mosaic_4326)  # returns degrees
lon_deg <- res_vals[1]
lat_deg <- res_vals[2]

#  set extent from a relevant sf layer
zoom_ext_4326 <- st_bbox(cypress_sf_4326)
bb <- zoom_ext_4326
bb["xmin"] <- bb["xmin"] - 0.2
bb["xmax"] <- bb["xmax"] + 0.2
bb["ymin"] <- bb["ymin"] - 0.2
bb["ymax"] <- bb["ymax"] + 0.2
xlim <- c(bb["xmin"], bb["xmax"])
ylim <- c(bb["ymin"], bb["ymax"])

# Plot
saplot_hu = ggplot() +
  # Basemap
  annotation_map_tile(type = "cartolight", zoomin = -1) +
  
  geom_tile(
    data = mosaic_df_sa,
    aes(x = x, y = y, fill = housing_units),
    width  = lon_deg * 1.001,
    height = lat_deg * 1.001,
    color  = NA,
    size   = 0,
    alpha  = 1.0
  ) +
  
  # Service area outline
  geom_sf(
    data  = cypress_sf_4326,
    aes(color = "Service Area"),
    fill      = NA,
    linewidth = 1
  ) +
  
  # Bank outline + fill by housing_units
  geom_sf(
    data = cypress_bank_filled,
    aes(
      fill  = housing_units,
      color = "Bank, avg. value 2021"
    ),
    linewidth = 0.4
  ) +
  
  # Viridis fill scale
  scale_fill_viridis_c(
    option   = "viridis",
    name     = "Downstream Housing Units",
    labels   = scales::comma_format(),
    direction = 1
  ) +
  
  # Manual color scale for outlines
  scale_color_manual(
    name   = NULL,
    values = c(
      "Bank, avg. value 2021" = "red",
      "Service Area"          = "steelblue"
    ),
    guide = guide_legend(
      override.aes = list(
        fill      = NA,
        linewidth = 1
      ),
      order = 2  # This sets the order - higher number means lower position
    )
  ) +
  
  # Coordinate limits
  coord_sf(
    xlim   = xlim,
    ylim   = ylim,
    expand = FALSE,
    crs    = st_crs(4326)
  ) +
  
  # Scale bar and north arrow
  annotation_scale(location = "bl", width_hint = 0.25) +
  annotation_north_arrow(
    location = "tr",
    which_north = "true",
    pad_x = unit(0.2, "in"),
    pad_y = unit(0.3, "in"),
    style = north_arrow_fancy_orienteering
  ) +
  
  # Theme + labels
  theme_minimal(base_size = 13) +
  labs(
    title    = "Remaining Wetlands (2021) in Big Cypress Mitigation Bank Service Area",
    subtitle = "Flood protection value of remaining 2021 wetlands, based on housing units downstream",
    fill     = "Downstream Housing Units",
    x        = NULL,
    y        = NULL
  ) +
  theme(
    legend.position = "right",
    plot.title      = element_text(size = 13, hjust = 0.5),
    plot.subtitle   = element_text(size = 12, hjust = 0.5),
    legend.title    = element_text(size = 11),
    legend.text     = element_text(size = 11)
  ) +
  guides(
    fill = guide_colorbar(order = 1)  # This sets the housing units legend order - lower number means higher position
  )

saplot_hu


########### PLOT SAME SA map without basemap #############

saplot_nbm = ggplot() +
  # Draw ocean *first* with soft blue
  annotate(
    "rect",
    xmin = xlim[1], xmax = xlim[2],
    ymin = ylim[1], ymax = ylim[2],
    fill = "aliceblue"
  ) +
  
  # Then draw land on top
  geom_sf(data = states, fill = "gray98", color = "gray70", linewidth = 0.2) +
  
  # Service area outline (legend label "Service Area")
  geom_sf(
    data  = cypress_sf_4326,
    aes(color = "Service Area"),
    fill = NA,
    linewidth = 1
  ) +
  
  # Bank outline (legend label "Bank, avg. value 2021") + fill by housing_units
  geom_sf(
    data = cypress_bank_filled,
    aes(
      fill  = housing_value,          # <— Map the same fill you want to see in the legend
      color = "Bank, avg. value 2021" # <— The outline color goes to a separate legend item
    ),
    linewidth = 0.4
  )+
  
  # Downstream housing units (raster layer)
  geom_tile(
    data   = mosaic_df_sa,
    aes(x = x, y = y, fill = housing_value),
    width  = lon_deg,
    height = lat_deg,
    alpha  = 0.8
  ) +
  
  # Only include the three specified cities
  geom_sf_label(
    data = places_fl %>% filter(NAME == "Cape Coral"),
    aes(label = NAME),
    stat = "sf_coordinates",
    nudge_y = 0.08,
    size = 3,
    label.size = 0.2,
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
    name = "Downstream Total Housing Value",
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
  
  coord_sf(xlim = xlim, ylim = ylim, expand = FALSE, crs = st_crs(target_crs)) +
  annotation_scale(location = "bl", width_hint = 0.25) +
  annotation_north_arrow(location = "tr", which_north = "true",
                         pad_x = unit(0.2, "in"), pad_y = unit(0.3, "in"),
                         style = north_arrow_fancy_orienteering) +
  
  theme_minimal(base_size = 13) +
  
  labs(
    title    = "Big Cypress Mitigation Bank Service Area",
    subtitle = "Remaining Wetlands 2021",
    fill     = "Downstream Total Housing Value", 
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

saplot_nbm


############### 

# Get downstream watersheds and centroids of loss and sa

###############


lossmosaic = rast("CONUS/Wetland Loss/Service Area Mosaics Loss - 3.5/Big_Cypress_MB_Phase_I-V_lossmosaic.tif")
samosaic = rast("CONUS/Service Area Mosaics 2021/Big_Cypress_MB_Phase_I-V_mosaic.tif")

get_downstream_hucs_from_mosaic <- function(mosaic_path, mosaic_id, sas, huc12_boundaries, htdh_path) {
  # Read the raster
  mosaic <- rast(mosaic_path)
  
  # Filter and reproject the service area polygon
  cypress <- sas[sas$ID == mosaic_id, ]
  cypress <- project(cypress, crs(mosaic))
  cypress_sf <- st_as_sf(cypress)
  
  # Reproject HUC12s to match raster
  huc12_boundaries <- st_transform(huc12_boundaries, st_crs(mosaic))
  
  # Intersect with HUC12s
  cypress_huc12s <- st_intersection(huc12_boundaries, st_union(cypress_sf))
  cypress_huc12s <- cypress_huc12s[, c("huc12", "geom")]
  
  # First band of raster
  loss_rast <- mosaic[[1]]
  
  # Convert to data.frame
  loss_vals <- as.data.frame(loss_rast, xy = TRUE, na.rm = TRUE)
  
  # Centroid
  mean_x <- mean(loss_vals$x)
  mean_y <- mean(loss_vals$y)
  
  loss_centroid <- st_as_sf(data.frame(x = mean_x, y = mean_y),
                            coords = c("x", "y"),
                            crs = crs(loss_rast))
  
  loss_centroid <- st_transform(loss_centroid, st_crs(cypress_huc12s))
  loss_centroid$label <- "Wetland Loss Centroid"
  
  # Find HUC12 for centroid
  loss_huc12 <- st_join(loss_centroid, cypress_huc12s, join = st_within)
  origin_huc12 <- loss_huc12$huc12
  
  # Load downstream mapping
  htdh <- read.csv(htdh_path)
  htdh$huc12 <- str_pad(htdh$huc12, width = 12, pad = "0")
  
  downstream_str <- htdh$downstream_hucs[htdh$huc12 == origin_huc12]
  downstream_ids <- strsplit(downstream_str, ",")[[1]]
  downstream_ids <- trimws(downstream_ids)
  downstream_ids <- downstream_ids[grepl("^\\d{12}$", downstream_ids)]
  
  downstream_hucs_sf <- cypress_huc12s %>% 
    filter(huc12 %in% downstream_ids)
  
  downstream_hucs_sf$source <- "Downstream Watersheds"
  
  return(list(
    loss_centroid = loss_centroid,
    downstream_hucs = downstream_hucs_sf
  ))
}

results_loss <- get_downstream_hucs_from_mosaic(
  mosaic_path = "CONUS/Wetland Loss/Service Area Mosaics Loss - 3.5/Big_Cypress_MB_Phase_I-V_lossmosaic.tif",
  mosaic_id = "Big_Cypress_MB_Phase_I-V",
  sas = sas,
  huc12_boundaries = readRDS("CONUS/Variables/huc12_boundaries.rds"),
  htdh_path = "CONUS/Variables/htdh.csv"
)

# Access outputs
loss_centroid <- results_loss$loss_centroid
loss_downstream_hucs <- results_los$downstream_hucs

#------------------------------ Plot just housing units demographic data --------------------------
library(tidycensus)
library(rnaturalearth)

options(tigris_use_cache = TRUE)
fl_hu_tracts<- get_acs(
  geography = "tract", 
  variables = "B25001_001",  # Number of housing units
  state = "FL",             
  year = 2021,               
  geometry = TRUE
) %>%
  rename(housing_units = estimate)
fl_hu_tracts = st_transform(fl_hu_tracts,st_crs(cypress_sf_4326))


cypress <- st_transform(cypress_sf, st_crs(cypress_sf_4326))
cypress_tracts <- fl_hu_tracts[st_intersects(fl_hu_tracts, cypress_sf_4326, sparse = FALSE), ]
# cypress_tracts <- st_intersection(fl_hu_tracts, cypress_sf_4326) #If we want to snip to SA boundary
cypress_tracts <- st_filter(fl_hu_tracts, cypress_sf_4326, .predicate = st_intersects)

# Get all census tracts

if (any(is.na(c(xlim, ylim)))) stop("xlim or ylim contains NA")

# Construct bbox safely
bbox_poly <- st_sfc(
  st_polygon(list(matrix(c(
    xlim[1], ylim[1],
    xlim[1], ylim[2],
    xlim[2], ylim[2],
    xlim[2], ylim[1],
    xlim[1], ylim[1]  # close the polygon
  ), ncol = 2, byrow = TRUE))),
  crs = st_crs(cypress_sf_4326)
)
map_tracts <- st_filter(fl_hu_tracts, bbox_poly, .predicate = st_intersects)

cities <- st_read("CONUS/Geo Boundaries/tl_2024_12_place/tl_2024_12_place.shp") %>%
  filter(NAME %in% c("Naples", "Fort Myers", "Cape Coral")) %>%
  mutate(geometry = st_centroid(geometry))  # keeps attributes + replaces geometry with centroid
cities <- st_transform(cities, st_crs(target_crs))

hu_tracts = ggplot() +
  annotation_map_tile(type = "cartolight", zoomin = -1) +
  # Tract-level housing units
  geom_sf(data = map_tracts, aes(fill = housing_units), color = NA, alpha = 0.9) +
  # Service area outline with color mapped to legend
  geom_sf(data = cypress_sf_4326, aes(color = "Service Area"), fill = NA, linewidth = 1) +
  
  # Bank outline LAST to bring it on top
  geom_sf(data = cypress_bank_4326, aes(color = "Bank"), fill = NA, linewidth = 1) +
  
  # Downstream HUCs - with aes() for legend
  geom_sf(data = downstream_hucs_sf, aes(color = "Downstream HUC12s"), fill = "black", alpha = 0.3, linewidth = 0.6) +
  
  # Loss centroid - with aes() for legend
  geom_sf(data = loss_centroid, aes(color = "Loss Centroid"), shape = 16, size = 3) +
  geom_sf(data = cities, aes(color = "City"), shape = 21, size = 2, fill = "white", stroke = 0.6) +
  geom_sf_text(data = cities, aes(label = NAME), size = 3, color = "darkgreen", fontface = "bold", nudge_y = 0.03)+
  # Fill scale for housing units
  scale_fill_distiller(
    palette = "YlGnBu",
    direction = 1,
    name = "Housing Units",
    labels = scales::comma_format()
  )+
  
  # Color scale for outlines - updated with new items
  scale_color_manual(
    name = NULL,
    values = c(
      "Bank" = "red",
      "Service Area" = "steelblue",
      "Downstream HUC12s" = "black",
      "Loss Centroid" = "black"
    )
  ) +
  
  # Legend order: fill first, then color
  guides(
    fill = guide_colorbar(order = 1),
    color = guide_legend(order = 2)
  ) +
  
  coord_sf(xlim = xlim, ylim = ylim, expand = FALSE, crs = st_crs(target_crs)) +
  # Scale bar and north arrow
  annotation_scale(location = "bl", width_hint = 0.25) +
  annotation_north_arrow(
    location = "tr",
    which_north = "true",
    pad_x = unit(0.2, "in"),
    pad_y = unit(0.3, "in"),
    style = north_arrow_fancy_orienteering
  ) +
  theme_minimal(base_size = 13) +
  
  labs(
    title    = "Big Cypress Mitigation Bank Service Area",
    subtitle = "Tract-level Housing Units (ACS 2021 5-Year Estimates)",
    fill     = "Housing Units",
    x        = NULL,
    y        = NULL
  ) +
  
  theme(
    legend.position   = "right",
    plot.title        = element_text(size = 13, hjust = 0.5),
    plot.subtitle     = element_text(size = 12, hjust = 0.5),
    legend.title      = element_text(size = 11),
    legend.text       = element_text(size = 11)
  )

hu_tracts

# Same plot but for remaining SA wetlands centroids --------------------

results_sa=get_downstream_hucs_from_mosaic(
  mosaic_path = "CONUS/Service Area Mosaics 2021/Big_Cypress_MB_Phase_I-V_mosaic.tif",
  mosaic_id = "Big_Cypress_MB_Phase_I-V",
  sas = sas,
  huc12_boundaries = readRDS("CONUS/Variables/huc12_boundaries.rds"),
  htdh_path = "CONUS/Variables/htdh.csv"
)
sa_centroid <- results_sa$loss_centroid
sa_downstream_hucs <- results_sa$downstream_hucs

hu_tracts_sa = ggplot() +
  annotation_map_tile(type = "cartolight", zoomin = -1) +
  # Tract-level housing units
  geom_sf(data = map_tracts, aes(fill = housing_units), color = NA, alpha = 0.9) +
  # Service area outline with color mapped to legend
  geom_sf(data = cypress_sf_4326, aes(color = "Service Area"), fill = NA, linewidth = 1) +
  
  # Bank outline LAST to bring it on top
  geom_sf(data = cypress_bank_4326, aes(color = "Bank"), fill = NA, linewidth = 1) +
  
  # Downstream HUCs - with aes() for legend
  geom_sf(data = sa_downstream_hucs, aes(color = "Downstream HUC12s"), fill = "black", alpha = 0.3, linewidth = 0.6) +
  
  # Loss centroid - with aes() for legend
  geom_sf(data = sa_centroid, aes(color = "Service Area Wetland Centroid"), shape = 16, size = 3) +
  geom_sf(data = cities, aes(color = "City"), shape = 21, size = 2, fill = "white", stroke = 0.6) +
  geom_sf_text(data = cities, aes(label = NAME), size = 3, color = "darkgreen", fontface = "bold", nudge_y = 0.03)+
  # Fill scale for housing units
  scale_fill_distiller(
    palette = "YlGnBu",
    direction = 1,
    name = "Housing Units",
    labels = scales::comma_format()
  )+
  
  # Color scale for outlines - updated with new items
  scale_color_manual(
    name = NULL,
    values = c(
      "Bank" = "red",
      "Service Area" = "steelblue",
      "Downstream HUC12s" = "black",
      "Service Area Wetland Centroid" = "black"
    )
  ) +
  
  # Legend order: fill first, then color
  guides(
    fill = guide_colorbar(order = 1),
    color = guide_legend(order = 2)
  ) +
  
  coord_sf(xlim = xlim, ylim = ylim, expand = FALSE, crs = st_crs(target_crs)) +
  # Scale bar and north arrow
  annotation_scale(location = "bl", width_hint = 0.25) +
  annotation_north_arrow(
    location = "tr",
    which_north = "true",
    pad_x = unit(0.2, "in"),
    pad_y = unit(0.3, "in"),
    style = north_arrow_fancy_orienteering
  ) +
  theme_minimal(base_size = 13) +
  
  labs(
    title    = "Big Cypress Mitigation Bank Service Area",
    subtitle = "Tract-level Housing Units (ACS 2021 5-Year Estimates)",
    fill     = "Housing Units",
    x        = NULL,
    y        = NULL
  ) +
  
  theme(
    legend.position   = "right",
    plot.title        = element_text(size = 13, hjust = 0.5),
    plot.subtitle     = element_text(size = 12, hjust = 0.5),
    legend.title      = element_text(size = 11),
    legend.text       = element_text(size = 11)
  )

hu_tracts_sa



# ECDF ---------------------------------

ecdf_loss = readRDS("CONUS\\Wetland Loss\\Extractions\\Extract Service Area Loss\\SA Summaries Loss\\SA ECDF Functions - Loss\\ecdf_fn_Big_Cypress_MB_Phase_I-V.rds")
ecdf_sa = readRDS("L:\\Wetland Flood Mitigation\\Extractions\\Extract Service Areas 2021\\SA ECDF Functions\\ecdf_fn_Big_Cypress_MB_Phase_I-V.rds")
ecdf_bank = readRDS("L:\\Wetland Flood Mitigation\\Extractions\\Extract Banks 2021\\Bank ECDF Functions\\bank_ecdf_fn_Big_Cypress_MB_Phase_I-V.rds")


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


ggplot(df_hu, aes(x = value, y = ecdf, color = source)) +
  geom_line(size = 1.5) +
  scale_x_continuous(labels = comma) +
  scale_color_discrete(
    labels = c(
      "Bank" = "Bank",
      "Service Area" = "Remaining Service Area (2021)",
      "Loss" = "Service Area Loss (1985–2021)"
    )
  ) +
  labs(
    title = "Downstream Housing Units of Wetlands – ECDF by Area",
    subtitle = "Big Cypress Wetland Mitigation Bank",
    x = "Downstream Housing Units", y = "Cumulative Probability",
    color=NULL
    
  ) +
  theme_minimal(base_size = 13) +
  theme(legend.position = "right",
        plot.title = element_text(hjust = 0.5),
        plot.subtitle = element_text(hjust = 0.5))


# HV ---

df_hv = df_all %>% 
  filter(variable=="Housing Value")

ggplot(df_hv, aes(x = value, y = ecdf, color = source)) +
  geom_line(size = 1.5) +
  scale_x_continuous(labels = comma) +
  scale_color_discrete(
    labels = c(
      "Bank" = "Bank",
      "Service Area" = "Remaining Service Area (2021)",
      "Loss" = "Service Area Loss (1985–2021)"
    )
  ) +
  labs(
    title = "Downstream Total Housing Value of Wetlands – ECDF by Area",
    subtitle = "Big Cypress Wetland Mitigation Bank",
    x = "Downstream Total Housing Value", y = "Cumulative Probability",
    color=NULL
    
  ) +
  theme_minimal(base_size = 13) +
  theme(legend.position = "right",
        plot.title = element_text(hjust = 0.5),
        plot.subtitle = element_text(hjust = 0.5))

# POP ---

df_pop = df_all %>% 
  filter(variable=="Population")

ggplot(df_pop, aes(x = value, y = ecdf, color = source)) +
  geom_line(size = 1.5) +
  scale_x_continuous(labels = comma) +
  scale_color_discrete(
    labels = c(
      "Bank" = "Bank",
      "Service Area" = "Remaining Service Area (2021)",
      "Loss" = "Service Area Loss (1985–2021)"
    )
  ) +
  labs(
    title = "Downstream Population of Wetlands – ECDF by Area",
    subtitle = "Big Cypress Wetland Mitigation Bank",
    x = "Downstream Population", y = "Cumulative Probability",
    color=NULL
    
  ) +
  theme_minimal(base_size = 13) +
  theme(legend.position = "right",
        plot.title = element_text(hjust = 0.5),
        plot.subtitle = element_text(hjust = 0.5))
########## Plot distribution ################

cypress_sa_extract = readRDS("Extractions/Extract Service Areas 2021/sa_extract_Big_Cypress_MB_Phase_I-V.rds")
cypress_bank_extract1 = readRDS("Extractions/Extract Banks 2021/bank_extract_Big_Cypress_MB_Phase_I-V.rds")
cypress_bank_extract2 = readRDS("Extractions/Extract Banks 2021/bank_extract_Big_Cypress_MB_Phase_VI.rds")
cypress_bank_extract = rbind(cypress_bank_extract1, cypress_bank_extract2)
cypress_loss_extract = readRDS("CONUS\\Wetland Loss\\Extractions\\Extract Service Area Loss\\sa_extract_Big_Cypress_MB_Phase_I-V.rds")

colnames(cypress_sa_extract)[2:4] <- c("housing_units", "housing_value", "population")
colnames(cypress_bank_extract)[2:4] <- c("housing_units", "housing_value", "population")
colnames(cypress_loss_extract)[2:4] <- c("housing_units", "housing_value", "population")

cypress_sa_extract$Source <- "Remaining Service Area (2021)"
cypress_bank_extract$Source <- "Bank"
cypress_loss_extract$Source <- "Service Area Loss (1985-2021)"

# Combine into one data frame
df <- bind_rows(cypress_sa_extract, cypress_bank_extract)
df <- bind_rows(df, cypress_loss_extract)

# Plot the distribution of log(housing_units) by source
ggplot(df, aes(x = log1p(housing_units), fill = Source, color = Source)) +
  geom_density(alpha = 0.3, size = 1.5) +
  scale_color_discrete(
    labels = c(
      "Bank" = "Bank",
      "Service Area" = "Remaining Service Area (2021)",
      "Loss" = "Service Area Loss (1985–2021)"
    )
  ) +
  scale_fill_discrete(
    labels = c(
      "Bank" = "Bank",
      "Service Area" = "Remaining Service Area (2021)",
      "Loss" = "Service Area Loss (1985–2021)"
    )
  ) +
  labs(
    title = "Distribution of Downstream Housing Units",
    subtitle = "Big Cypress Wetland Mitigation Bank",
    x = "Log Downstream Housing Units", y = "Density",
    color = NULL,
    fill = NULL
  ) +
  theme_minimal(base_size = 13) +
  theme(
    plot.title = element_text(hjust = 0.5),
    plot.subtitle = element_text(hjust = 0.5),
    legend.position = "right"
  )

# HV 

ggplot(df, aes(x = log1p(housing_value), fill = Source, color = Source)) +
  geom_density(alpha = 0.3, size = 1.5) +
  scale_color_discrete(
    labels = c(
      "Bank" = "Bank",
      "Service Area" = "Remaining Service Area (2021)",
      "Loss" = "Service Area Loss (1985–2021)"
    )
  ) +
  scale_fill_discrete(
    labels = c(
      "Bank" = "Bank",
      "Service Area" = "Remaining Service Area (2021)",
      "Loss" = "Service Area Loss (1985–2021)"
    )
  ) +
  labs(
    title = "Distribution of Downstream Total Housing Value",
    subtitle = "Big Cypress Wetland Mitigation Bank",
    x = "Log Downstream Total Housing Value", y = "Density",
    color = NULL,
    fill = NULL
  ) +
  theme_minimal(base_size = 13) +
  theme(
    plot.title = element_text(hjust = 0.5),
    plot.subtitle = element_text(hjust = 0.5),
    legend.position = "right"
  )

# POP 
ggplot(df, aes(x = log1p(population), fill = Source, color = Source)) +
  geom_density(alpha = 0.3, size = 1.5) +
  scale_color_discrete(
    labels = c(
      "Bank" = "Bank",
      "Service Area" = "Remaining Service Area (2021)",
      "Loss" = "Service Area Loss (1985–2021)"
    )
  ) +
  scale_fill_discrete(
    labels = c(
      "Bank" = "Bank",
      "Service Area" = "Remaining Service Area (2021)",
      "Loss" = "Service Area Loss (1985–2021)"
    )
  ) +
  labs(
    title = "Distribution of Downstream Population",
    subtitle = "Big Cypress Wetland Mitigation Bank",
    x = "Log Downstream Population", y = "Density",
    color = NULL,
    fill = NULL
  ) +
  theme_minimal(base_size = 13) +
  theme(
    plot.title = element_text(hjust = 0.5),
    plot.subtitle = element_text(hjust = 0.5),
    legend.position = "right"
  )

################# T-tests ######################

library(dplyr)

cypress_subset <- df %>%
  filter(Source %in% c("Bank", "Service Area Loss (1985-2021)")) %>%
  mutate(log_pop = log1p(population))

t.test(population ~ Source, data = cypress_subset)
t.test(log_pop ~ Source, data = cypress_subset)
