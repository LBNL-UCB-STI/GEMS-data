# CODE TO COMPILE AND CLEAN ALL INPUTS FOR SECOND STAGE CLUSTERING 
# UPDATED SECOND STAGE FOR THE TRANSP GEO RE-SUBMISSION
# need to make one observation per spatial_id (cbsa or county)
# NOV 30 2020: added urban divisions
# NATALIE POPOVICH
################
# ADDED temprorary environment -- can be removed once all the R packages are installed
# mywd = "C:/FHWA/For FHWA folks/TransportGeography-Revision"
# setwd(mywd)
# 
# figuredir <- "./Figures"
# inputsdir <- "./Data/InputData"
# datadir <- "./Data/OutputData"
# rdatadir <-  "./Data/RData"
# tabdir <- './Tables'
# resultsdir <- './Results'
# 
# 
# library(sf)
# library(dplyr)
# library(stargazer)
# library(psych)
# library(data.table)
# library(tigris)
# library(tidyr)
# library(openxlsx)

# Identifiers and urban population divisions and IDs
div <- read.csv(file = file.path(inputsdir,"urban_divisions.csv")) %>%
  select(spatial_id, GEOID, fhwa_type, fhwa_type_num) %>%
  rename(tract = GEOID)

#########################
# INPUTS REQUIRED:
# 1: percent of each microtype in each spatial ID (CBSA or county)
# 2: HHI of commutes. distribution of PMT throughout census tracts
# 3: number of employment centers (polycentricity metric)
# 4: urban divisions
##################

# 1: microtype distribution in each spatial ID (Ling's new clusters Nov 20 2020)
micro <- read.csv(file = file.path(inputsdir,"clustering_outputs.csv"))  %>%
  select(tract, cluster6) # 6 total microtypes

dat <- micro %>% 
  right_join(div, by = "tract") %>%
  group_by(spatial_id, cluster6) %>% # generate percentage of tracts of each type
  summarise(n = n()) %>%
  mutate(pct = n/sum(n)) %>%
  select(-n)

#reshape wide to turn into individual inputs, should have 2250 total rows
wide <- dat %>%
  spread(cluster6, pct)

# replace NAs with zeros
wide[is.na(wide)] <- 0

wide <- wide %>% 
  rename(na = "<NA>") %>% 
  rename_at(vars(2:8) ,function(x){paste0("pct_micro.", x)}) %>% ungroup()

########
#2) commute metrics: HHI
##############
hhi <-read.csv(file = file.path(inputsdir,"regional_metrics_v2.csv")) %>%
  select(spatial_id, hhi_normalized)

dat <- wide %>% 
  merge(hhi, by = "spatial_id") 

#############
# 3) EMPLOYMENT POLYCENTRICITY
# Note this was done in QGIS, but could try to do it in R
##############
# Steps:
# 1 create one shapefile by spatial_id: county and CBSAs in one layer
# 2 read in shapefile of microtypes
# 3 dissolve polygons with same microtype cluster number that are adjacent to each other
# 4 assign centroid to each polygon part
# 5 count # points in each spatial ID (cbsa or county)

# Step 1: create one shapefile by spatial_id: county and CBSAs in one layer
# read in boundary shapefiles for region clusters
library(sf)
cty <- st_read(file.path(inputsdir, "Counties/combined_counties.shp")) %>% 
  select(GEOID, CBSAFP, ALAND, AWATER) %>% 
  mutate(spatial_id = case_when(
  is.na(CBSAFP) ~ GEOID,
  TRUE ~ CBSAFP)) 

# exporting shapefile to dissolve by spatial ID in QGIS
st_write(cty, file.path(datadir, "spatial_ids.shp"), driver = "Esri Shapefile", append=FALSE)

#[DISSOLVE DONE IN QGIS BY FIELD = spatial_id]
#bring in combined shapefile with polygons attributes by spatial ID
regions <- st_read(file.path(inputsdir, "spatial_ids_dissolved/spatial_ids_dissolved.shp"))

# Step 2 read in shapefile of microtypes (probably easiest in QGIS)
tracts <- st_read(file.path(inputsdir, "combined_tracts/combined_tracts.shp")) %>%
  mutate(tract = as.numeric(GEOID))

# merge with cluster labels 
tracts <- tracts %>% 
  merge(micro, by = "tract")

st_write(tracts, file.path(datadir,"tract_clusters.shp"), driver = "Esri Shapefile", append=FALSE)

#############################
# Step 3-5  [DONE IN QGIS]
# dissolve polygons with same microtype cluster number that are adjacent 
                #to each other ("Dissolve by field = cluster6")
# assign centroid to each polygon part ("Centroids")
# filter by cluster6 only ("Extract by Attribute")
# count # points in each spatial ID (cbsa or county) ("Count Points in Polygons")
            # layers: tract_clusters.shp and spatial_ids_dissolved.shp
# export csv with number of points in each spatial id
#################################

# read results back in
count <- read.csv(file = file.path(inputsdir,"employment_centers_transp_geo.csv")) %>%
  select(spatial_id, NUMPOINTS) %>%
  rename(emp_centers = NUMPOINTS)

dat <- dat %>%
  merge(count, by = "spatial_id", all.x = T)

# should have a total of 1,435 employment centers (results of merged cluster6 = 6)
sum(dat$emp_centers) # we do

# reading in urban boundaries, need to assign each spatial region to one category
# for counties/cbsas with tracts belonging to multiple definitions, round up
bounds <- div %>% 
  select(spatial_id, fhwa_type, fhwa_type_num) %>% 
  distinct() %>%
  group_by(spatial_id) %>%
  mutate(max = max(fhwa_type_num)) %>%
  filter(fhwa_type_num == max | is.na(fhwa_type_num)) %>%
  select(-max)

dat <- dat %>%
  merge(bounds, by = "spatial_id")

##############
# DATA CLEANING
#################
colSums(is.na(dat)) # check number of missing values for each variables

write.csv(dat, file = file.path(inputsdir, "geotypes_inputs_transp_geo_060623.csv"), row.names = F)
