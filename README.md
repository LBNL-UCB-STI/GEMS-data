# GEMS Data Generation 

## POC: Xiaodan Xu, Ph.D. (XiaodanXu@lbl.gov)
## Latest update: 12/07/2023


## Theme A: Geographic boundary
### a1. collecting micro-geotype boundary
step 1: collecting and cleaning spatial crosswalk file for LODES8 (using 2020 census boundary)
**code**: [geography/clean_lodes8_crosswalk.py](clean_lodes8_crosswalk.py)
**input**: downloaded state-by-state crosswalk from LEHD website (the R API doesn't work)
https://lehd.ces.census.gov/data/
**output**:spatial_boundary/CleanData/cleaned_lodes8_crosswalk.csv

step 2: Collecting Census 2020 Tiger line/boundary for census tract, county and CBSA
**code**:[geography/0_clean_boundaries.R](0_clean_boundaries.R)
**Input**: queried 2020 Census boundary from API in R;
Cleaned census crosswalk file downloaded here:
https://www2.census.gov/geo/docs/maps-data/data/rel2020/tract/tab20_tract20_tract10_natl.txt
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
**output**:
* Demography/CleanData/acs_data_tracts_{date}.csv


## Theme C: Collect demand related attributes
### c1. clean and processing LEHD LODES8 data
step 1: collecting LEHD LODES 8 data
**code**:[0_clean_lehd_2017.R](demand/0_clean_lehd_2017.R)
**Input**: query latest LEHD LODES8 data from API for all states at tract level 
(2021 for most states, with AK, AR, MS only has older data available)
* load and clean LODES8 Workplace Area Characteristics (WAC) at census tract level
* load and clean LODES8 Origin-Destination (OD) data at census tract level, include auxilliary data for cross-state commute
**output**:
* Demand/CleanData/wac_tract_{year}.csv
* Demand/CleanData/OD/*

step 2: calculate commute distance (using great circle distance between OD centroids)
**code**:[generate_od_distance.py](demand/generate_od_distance.py)
**Input**: 
Demand/CleanData/OD/* and spatial_boundary/CleanData/combined_tracts_{year}.csv
(2021 for most states, with AK, AR, MS only has older data available)
**output**:
* Demand/CleanData/OD_distance/*

## Theme D: Network generation
### d1. processing OSMNX data at census tract level
**code**: [generate_OSMNX_metrics.py](network/generate_OSMNX_metrics.py)
**Input**: queried 2023 OSM metrics from API in R and spatial boundary from step a1 above
* load network statistics from OSMNX
* if no statistics found, fill in tract-level attributes with NA
**output**:
* Network/CleanData/OSMNX/*

## Theme X: spatial clustering
### X1. compile demand attributes at census tract level
**code**: [demand_variable_generation.py](spatial_cluster/demand_variable_generation.py)
**Input**: compile inputs from various themes, including:
* load spatial boundary from Theme A: spatial_boundary/CleanData/combined_tracts_{year}.geojson and csv
* load population data from Theme B: Demography/CleanData/acs_data_tracts_{date}.csv
* load lehd wac data from Theme C: Demand/CleanData/wac_tract_{year}.csv AND Demand/CleanData/OD_distance/*
[Placeholder for more attributes]
**output**:
* Demand/CleanData/microtype_inputs_demand.csv