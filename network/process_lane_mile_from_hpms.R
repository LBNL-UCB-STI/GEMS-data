library(sf)
library(dplyr)
library(data.table)
library(tidycensus)
library(tools)

readRenviron("~/.Renviron")

merge_census_hpms<- function(state='california2017',state_tracts,year){
  print(state)
  hpms_geometry <- st_read(paste0(state,'/',toTitleCase(substr(f,10,50)), '.shp')) # read HPMS data
  crs_hpms <- st_crs(hpms_geometry)
  
  # assign tracts to each link
  hpms_geometry_by_tracts = st_intersection(st_zm(hpms_geometry), state_tracts)
  
  hpms_geometry_by_tracts$Length = st_length(hpms_geometry_by_tracts) # re-generate link length after split network by tracts
  hpms_geometry_by_tracts <- hpms_geometry_by_tracts %>% mutate(lanemiles = as.numeric(Through_La * Length / 1609.34)) # compute lane miles
  hpms_geometry_by_tracts <- hpms_geometry_by_tracts %>% filter(lanemiles > 0) # remove invalid links
  hpms_geometry_by_tracts <- hpms_geometry_by_tracts %>% mutate(g_type = st_geometry_type(.)) 
  hpms_geometry_by_tracts <- hpms_geometry_by_tracts %>% filter(g_type %in% c('LINESTRING', 'MULTILINESTRING'))
  st_write(hpms_geometry_by_tracts, paste0( state, '_HPMS_with_',year,'_GEOID_LANEMILE.geojson'),append=FALSE) # save merged data
}

# load HPMS data
setwd("C:/FHWA_R2/Network/RawData/") # set input directory
file_list <- list.files("HPMS2017", full.names = TRUE) # list all HPMS data files
file_list <- file_list[!grepl("zip", file_list, ignore.case = TRUE)]

# load state census tracts (whole U.S.)
tracts <- st_read("C:/FHWA_R2/spatial_boundary/CleanData/combined_tracts_2020.geojson")
tracts <- st_transform(tracts, crs = crs_hpms)
sf::sf_use_s2(FALSE) 

for (f in file_list){
  merge_census_hpms(f,tracts,2020)
}


