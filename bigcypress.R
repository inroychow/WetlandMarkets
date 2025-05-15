
###################################################################################################################
######                       Case Study: Big Cypress MB in Southwest Florida                       ################
######     has the most lost wetland cells and the highest value lost in SA compared to bank.      ################
###################################################################################################################


setwd("C:/Users/indumati/Box/Paper2_final")

purrr::walk(c("future", "future.apply", "ggplot2", "tidyverse", "tibble", "progressr", 
              "data.table", "dplyr", "terra", "tigris", "sf", "scales", "ggspatial", 
              "rosm", "raster", "rnaturalearth"), library, character.only = TRUE)


loss_summary = readRDS("Wetland Loss/Extractions and Summaries/full_summary_loss.rds")
bank_summary = readRDS("Extractions and Summaries/Summaries Banks/full_summary_bank.rds")

sas= readRDS("Service Areas/ServiceAreas_agg.rds")
banks = readRDS("Bank Footprints/footprints_and_buffers.rds")


#Find which states intersect this SA ----------------------------------------

cypress_banks = banks[banks$Name %in% c("Big_Cypress_MB_Phase_I-V", "Big_Cypress_MB_Phase_VI")] # This bank has two parallel footprints
cypress_banks = terra::union(cypress_banks)
cypress_sa = sas[sas$ID == "Big_Cypress_MB_Phase_I-V",]

states <- ne_states(country = "United States of America", returnclass = "sf")
states_albers <- st_transform(states, crs(cypress_sa))
states_vect <- vect(states_albers)
intersected <- terra::intersect(states_vect, cypress_sa)
intersected$name #Florida

fl_sf <- subset(states_albers, name %in% "Florida")


# Load loss raster
cypress_mosaic_loss = rast("Wetland Loss/Service Area Mosaics Loss/Big_Cypress_MB_Phase_I-V_lossmosaic.tif")
mosaic_df_loss <- as.data.frame(cypress_mosaic_loss, xy = TRUE, na.rm = TRUE)
colnames(mosaic_df_loss)[3:5] <- c("housing_units", "housing_value", "population")

cypress_sa <- st_as_sf(cypress_sa)
cypress_banks = st_as_sf(cypress_banks)

# Note: ggspatial assumes map is in lon/lat (EPSG:4326) WGS8, so reproject everything
target_crs <- 4326

cypress_sa_4326 <- st_transform(cypress_sa, target_crs)
cypress_bank_4326 <- st_transform(cypress_banks, target_crs)
fl_sf_4326 <- st_transform(fl_sf, target_crs)
coords_4326 <- terra::project(cbind(mosaic_df_loss$x, mosaic_df_loss$y), crs(cypress_mosaic_loss), "EPSG:4326")

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
  raster_path = "Service Area Mosaics/Big_Cypress_MB_Phase_I-V_mosaic.tif",
  polygon_sf = cypress_bank_4326,
  variable = "housing_units"  # can also be "housing_value" or "population"
)


# ---------------------- Plot Loss WITHOUT basemap ----------------------

cell_res <- res(mosaic_raster_4326)  # returns a vector: c(xres, yres)
lon_deg <- cell_res[1] # For use in geom_tile
lat_deg <- cell_res[2]

places_fl <- places(state = "FL", year = 2021)  

lossplot_nbm = ggplot() +
  # Draw ocean with soft blue
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
      fill  = housing_units,          
      color = "Bank, avg. value 2021" 
    ),
    linewidth = 0.4
  )+
  
  # Downstream housing units (raster layer)
  geom_tile(
    data   = mosaic_df_4326_loss,  # Use the reprojected data
    aes(x = x, y = y, fill = housing_units),
    width  = lat_deg,  
    height = lon_deg,  
    alpha  = 0.8
  )+
  
  # Cape Coral label
  geom_sf_label(
    data = places_fl %>% filter(NAME == "Cape Coral"),
    aes(label = NAME),
    stat = "sf_coordinates",
    nudge_y = 0.08,
    size = 3,
    label.size = 0.3,
    label.padding = unit(1, "mm")
  ) +
  
  # Naples and Fort Myers labels
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
    # Override the default legend key order
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

#  -------------------------------- Plot remaining SAs mosaic -----------------------------------------------------

# Load the raster
cypress_mosaic_sa <- rast("Service Area Mosaics/Big_Cypress_MB_Phase_I-V_mosaic.tif")

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


# ECDF ---------------------------------

ecdf_loss = readRDS("ECDF Functions/Loss ECDF Functions/ecdf_fn_Big_Cypress_MB_Phase_I-V.rds")
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

cypress_sa_extract = readRDS("Extractions and Summaries/Extract Service Areas/sa_extract_Big_Cypress_MB_Phase_I-V.rds")
cypress_bank_extract1 = readRDS("Extractions and Summaries/Extract Banks/bank_extract_Big_Cypress_MB_Phase_I-V.rds")
cypress_bank_extract2 = readRDS("Extractions and Summaries/Extract Banks/bank_extract_Big_Cypress_MB_Phase_VI.rds")
cypress_bank_extract = rbind(cypress_bank_extract1, cypress_bank_extract2)
cypress_loss_extract = readRDS("Wetland Loss/Extractions and Summaries/Extract Service Area Loss/sa_extract_Big_Cypress_MB_Phase_I-V.rds")

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
