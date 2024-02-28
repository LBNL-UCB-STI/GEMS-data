# GEMS Data Generation 

## POC: Xiaodan Xu, Ph.D. (XiaodanXu@lbl.gov)
## Latest update: 02/27/2024


## Theme A: Geographic boundary
### a1. collecting micro-geotype boundary

**step 1: collecting and cleaning spatial crosswalk file for LODES8** (using 2020 census boundary)

**code**: [clean_lodes8_crosswalk.py](geography/clean_lodes8_crosswalk.py)

**input**: downloaded state-by-state crosswalk from LEHD website (the R API doesn't work)
https://lehd.ces.census.gov/data/

**process**: Clean up raw crosswalk files by state and create single lookup table for all states

**output**:spatial_boundary/CleanData/cleaned_lodes8_crosswalk.csv


**step 2: Collecting Census 2020 Tiger line/boundary for census tract, county and CBSA**

**code**:[0_clean_boundaries.R](geography/0_clean_boundaries.R)

**Input**: queried 2020 Census boundary from API in R;
Cleaned census crosswalk file downloaded here:
https://www2.census.gov/geo/docs/maps-data/data/rel2020/tract/tab20_tract20_tract10_natl.txt

**process**:
* load census tract geometry and land/water area (for microtypes)
* load and combine county and CBSA boundary
* load and clean spatial crosswalks (2020-2010 crosswalk, county-cbsa-tract crosswalk)

**output**:
* spatial_boundary/CleanData/combined_tracts_{year}.geojson and csv
* spatial_boundary/CleanData/combined_county_{year}.geojson and csv
* spatial_boundary/CleanData/combined_geotype_unit_{year}.geojson and csv
* spatial_boundary/CleanData/cleaned_lodes8_crosswalk_with_ID.csv

## Theme B: Demograhic characteristics
### b1. Collecting ACS data at census tract level

**code**: [1_ACS_compile_tracts.R](demographic/1_ACS_compile_tracts.R)

**Input**: query 2021 ACS 5-Year estimates from API

**process**: collecting tract-level demographic variables for persons, households and housing units from ACS 5-year estimates using tidycensus (with ACS API) 

**output**:
* Demography/CleanData/acs_data_tracts_{date}.csv


## Theme C: Collect demand related attributes
### c1. clean and processing LEHD LODES8 data

**step 1: collecting LEHD LODES 8 data**

**code**:[0_clean_lehd_2017.R](demand/0_clean_lehd_2017.R)

**Input**: query latest LEHD LODES8 data from API for all states at tract level 
(2021 for most states, with AK, AR, MS only has older data available)

**process**:
* load and clean LODES8 Workplace Area Characteristics (WAC) at census tract level using LEHDR (R package for accessing LEHD data)
* load and clean LODES8 Origin-Destination (OD) data at census tract level using LEHDR, include auxiliary data for cross-state commute

**output**:
* Demand/CleanData/wac_tract_{year}.csv
* Demand/CleanData/OD/*


**step 2: calculate commute distance** (using great circle distance between OD centroids)

**code**:[generate_od_distance.py](demand/generate_od_distance.py)

**Input**: 
Demand/CleanData/OD/* and spatial_boundary/CleanData/combined_tracts_{year}.csv
(2021 for most states, with AK, AR, MS only has older data available)

**process**: 
* Calculate euclidean distance between each OD pair
* For intrazonal OD, using 1/3* sqrt(area) as a proxy for distance

**output**:
* Demand/CleanData/OD_distance/*

## Theme D: Collect land use attributes
### D1. land use characteristics from NLCD data

**step 1: collecting and cleaning NLCD data**
**Code**:[process_NLCD_data.R](geography/process_NLCD_data.R)

**Input**:
* CONUS: Land_use/RawData/US_2020_nlcd_shapefiles_24mar2023/NLCD
* Alaska: Land_use/RawData/emiss_shp2017/NLCD
* Hawaii: Land_use/RawData/hi_hawaii_2010_ccap_hr_land_cover20150120

**Processes**:
* For CONUS and AK, assign centroids of grid cell to each census tract, and aggregate areas by land use types in each census tract
* For HI, generate grid cell (polygons) from pixels in raster image (island by island), and follow the similar process as CONUS and AK

**Output**:
* Land_use/CleanData/tract_level_land_use_no_ak.csv
* Land_use/CleanData/tract_level_land_use_ak.csv
* Land_use/CleanData/HI_NLCD_2010_{island_names}.csv

**step 2: Compile and impute NLCD data**

**Code**:[compile_nlcd_data.py](geography/compile_nlcd_data.py)

**Input**: 
* Land_use/CleanData/tract_level_land_use_no_ak.csv
* Land_use/CleanData/tract_level_land_use_ak.csv
* Land_use/CleanData/HI_NLCD_2010_{island_names}.csv
* spatial_boundary/CleanData/combined_tracts_{year}.geojson

**Processes**:
* Unify land use type names, e.g., map different types of imperious developed land to variable 'imperious developed'
* Combine land use data for all states
* Calculate fractions of imperious land and developed open space by census tract, and impute missing value using values from nearest census tracts

**Output**:
Land_use/CleanData/processed_NLCD_data.csv (for all land use types, such as forest, agriculture, developed)
Land_use/CleanData/imputed_NLCD_data_dev_only.csv


### D2. uban area definition from U.S.Census Bureau
**Code**:[1_clean_urban_areas_definitions.R](geography/1_clean_urban_areas_definitions.R)

**Input**: 
* Spatial crosswalk between urban area and census block: spatial_boundary/RawData/2020_UA_BLOCKS.txt'

**Processes**:
* Merge census urban area (UA) boundary with all census tract
* Generate urban/rural indicators based on merge results (1 - if within UA, 0 - otherwise)
* Assign urban area definition used in V1 typology for validation (majority of the classification results are the same)

**Output**:
spatial_boundary/CleanData/urban_divisions_2021.csv

## Theme E: Network generation
### e1. processing OSMNX data at census tract level

**code**: [generate_OSMNX_metrics.py](network/generate_OSMNX_metrics.py)

**Input**: queried 2023 OSM metrics from API in R and spatial boundary from step a1 above

* load network statistics from OSMNX
* if no statistics found, fill in tract-level attributes with NA

**output**:
* Network/CleanData/OSMNX/*

## Theme X: Spatial clustering
### X1. Develop socio-economic microtype

**step 1: compile demand attributes at census tract level**

**code**: [demand_variable_generation.py](spatial_cluster/demand_variable_generation.py)

**Input**: compile inputs from various themes, including:
* load spatial boundary from Theme A: spatial_boundary/CleanData/combined_tracts_{year}.geojson and csv
* load population data from Theme B: Demography/CleanData/acs_data_tracts_{date}.csv
* load lehd wac data from Theme C: Demand/CleanData/wac_tract_{year}.csv AND Demand/CleanData/OD_distance/*
* Load NLCD land use data from Theme D: Land_use/CleanData/imputed_NLCD_data_dev_only.csv
* Load urban area definition from Theme D: spatial_boundary/CleanData/urban_divisions_2021.csv

**Processes**:
* Load all data sources needed for socio-economic clusters
* Filter out census tracts with only water surface (no land)
* Generate tract-level variables using 2020 census boundary and count number of missing values for each variable generated

**output**:
* Demand/CleanData/microtype_inputs_demand.csv

**step 2: Develop and validate socio-economic microtype**

**code**: [demand_microtype_cluster.R](spatial_cluster/demand_microtype_cluster.R)

**dependency**:[initialization.R](spatial_cluster/initialization.R) and [functions.R](spatial_cluster/functions.R)

**Input**:
* Demand/CleanData/microtype_inputs_demand.csv
* spatial_boundary/CleanData/cleaned_lodes8_crosswalk_with_ID.csv

**Processes**:
* load data, performing cleaning, scaling, imputation and then split by rural/urban boundary
* spatial cluster for rural and urban respectively
* validate spatial clusters through visualization

**Output**:
* Demand/Results/microtypes_inputs_demand_scaled.csv
* Demand/Results/clustering_outputs_with_raw_data.csv