# Description of the data and file structure

We created a nationwide 30 m raster dataset to analyze flood protection value in wetland mitigation banks across the continental United States.

## Files and variables

### Data

- **footprints_and_buffers.shp** (or **footprints_and_buffers.rds**): Shapefile of mitigation bank footprint polygons.
- **ServiceAreas_agg.shp** (or **ServiceAreas_agg.rds**): Shapefile of mitigation bank service area polygons.
- **service_area_huc12_to_downstream_tracts.csv**: Raw output of `Upstream_iteration.ipynb`; one column contains HUC12 watersheds and the other contains a list of downstream census tracts.
- **htdt_conus.rds**: Dataframe relating HUC12 watershed boundaries to all respective downstream census tracts.
- **huc12_boundaries.gpkg**: Shapefile of HUC12 watershed polygons for the entire U.S.
- **LCMAP_CU_1985_V13_LCPRI.tif**: Raster of land cover data from the USGS LCMAP product for 1985.
- **LCMAP_CU_2021_V13_LCPRI.tif**: Raster of land cover data from the USGS LCMAP product for 2021.
- **places_fmv_vacant.tif**: Raster of vacant land value data from Nolte (2020).

### Code/software

We use **RStudio version 4.5.1** to run all code, with the exception of one script run in **Python version 3.13.0**.

## Script files

- **create_variables.R**  
  Creates and cleans necessary variables for analysis, including wetland raster generation, census data retrieval, and watershed-to-tract mapping.

- **Upstream_iteration.ipynb**  
  Maps HUC12 watersheds to their downstream census tracts.

- **run_conus.R**  
  Defines and applies a distance-decay function to wetland pixels across all mitigation banks, saving output as classified HUC12 watershed rasters.

- **missing_huc12s_analysis.R**  
  Applies the distance-decay function to missing HUC12s where downstream tract centroids were out of range, assuming uniform population distribution.

- **run_wetlandloss.R**  
  Identical to `run_conus.R` but applied to a raster of wetland loss (wetland-to-development conversion between 1985 and 2021).

- **mosaic.R**  
  Stitches together individual HUC12 classified rasters into larger mitigation bank service area rasters.

- **extract_cell_values.R**  
  Extracts classified cell values into dataframes for each bank footprint and service area. Also applies to the classified wetland loss raster.

- **summarize_ecdf_fn.R**  
  Computes empirical cumulative distribution functions (ECDFs) for extracted cell values in each service area, for both existing and lost wetlands.

- **pooling.R**  
  Creates pooling and pool mapping for mitigation banks with identical service area geometries.

- **generate_main_figures.R**  
  Generates the primary figures used in the manuscript.

- **generate_supp_figures.R**  
  Produces supplemental figures for the manuscript.
