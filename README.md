# GEMS-data

## Geographic boundary
### collecting micro-geotype boundary
step 1: collecting and cleaning spatial crosswalk file
code: [geography/clean_lodes8_crosswalk.py](clean_lodes8_crosswalk.py)

step 2: Collecting Census 2020 Tiger line/boundary for census tract, county and CBSA
code:[geography/0_clean_boundaries.R](0_clean_boundaries.R)

## Network generation
### processing OSMNX data at census tract level
code: [generate_OSMNX_metrics.py](network/generate_OSMNX_metrics.py)
* load network statistics from OSMNX
* if no statistics found, fillin tract-level attributes with NA