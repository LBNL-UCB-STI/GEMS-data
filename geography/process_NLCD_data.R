# NLCD data processing
# Created on 12/19/23 by Xiaodan Xu, create fraction of developed land and agriculture land
######################################################################
# LJ add: directories so that the files are saved in the right places
mywd <- "C:/FHWA_R2/"
setwd(mywd)

datadir <- "RawData"
cleandir <- "CleanData"

library(rgdal)
require(RCurl)
library(purrr) # LJ add: to use reduce
require(parallel)
library(sf)
library(sp)
library(rgeos)
library(maptools)
library(maps)
library(tigris)
library(stringr)
library(lwgeom) # to calculate area
library(tidyr) # to use spread func

sf_use_s2(FALSE)

# load NLCD data
nlcd_path <- 'US_2020_nlcd_shapefiles_24mar2023/NLCD/nlcd_2019_land_cover_l48_20210604_500m_ll.shp'
nlcd_data <- st_read(file.path('Land_use', datadir, nlcd_path))

# load census boundary
tract_path <- 'combined_tracts_2020.geojson'
us_census_tract <- st_read(file.path('spatial_boundary', cleandir, tract_path))
us_census_tract <- st_transform(us_census_tract, 4326) # change to NLCD CRS
print(st_crs(us_census_tract))

nlcd_data$area <- st_area(nlcd_data)
nlcd_data_centroid <- st_centroid(nlcd_data)

# nlcd_data_sample <- head(nlcd_data, 1000)
# plot(nlcd_data_sample[, 'GRIDCODE'])
 

# plot(nlcd_data_centroid_sample[, 'GRIDCODE'])

nlcd_data_centroid <- st_join(nlcd_data_centroid, us_census_tract,
                             join = st_nearest_feature, left = TRUE)

nlcd_data_centroid_df <- nlcd_data_centroid %>% st_drop_geometry()
write.csv(nlcd_data_centroid_df, file.path('Land_use', cleandir, 'assigned_nlcd_to_tract.csv'))
          
nlcd_data_centroid_sample <- head(nlcd_data_centroid, 1000)          
selected_tracts <- unique(nlcd_data_centroid_sample$GEOID)
selected_tracts_sf <- us_census_tract[us_census_tract$GEOID %in% selected_tracts,]
plot(st_geometry(selected_tracts_sf))
plot(nlcd_data_centroid_sample[, 'GEOID'], add= TRUE)

# nlcd_data_by_tract <- nlcd_data_centroid_df %>% 
#   dplyr::group_by(GEOID, GRIDCODE) %>% 
#   dplyr::summarise(area = sum(area)) %>% 
#   tidyr::spread(GRIDCODE, area)

nlcd_data_centroid_sample_df <- nlcd_data_centroid_sample %>% st_drop_geometry()

nlcd_data_centroid_df <- nlcd_data_centroid_df %>% 
  dplyr::mutate(area = as.numeric(area)) # unit in m^2

nlcd_data_by_tract<- nlcd_data_centroid_df %>%
  tidyr::pivot_wider(names_from = GRIDCODE, values_from = area, 
                   values_fn = sum, values_fill = 0)

write.csv(nlcd_data_by_tract, file.path('Land_use', cleandir, 'tract_level_land_use_no_ak.csv'))

# fill land use for AK with old NLCD
ak_hi_nlcd_path <- 'emiss_shp2017/NLCD/CONUS_AK_NLCD_2011_500m_WGS.shp'
ak_hi_nlcd_data <- st_read(file.path('Land_use', datadir, ak_hi_nlcd_path))

ak_hi_nlcd_data$area <- st_area(ak_hi_nlcd_data)
ak_nlcd_data_centroid <- st_centroid(ak_hi_nlcd_data)
ak_nlcd_data_sample <- head(ak_nlcd_data_centroid, 1000)
plot(ak_nlcd_data_sample[, 'GRIDCODE'])
# get ak boundary
states <- tigris::states(cb = TRUE)
ak_hi_boundary <- states %>% dplyr::filter(STUSPS == 'AK' | STUSPS == 'HI')
plot(st_geometry(ak_hi_boundary))
ak_hi_boundary <- st_transform(ak_hi_boundary, 4326) # re-project data

# using 'st_intersects', which gives a boolean for if points within polygon
ak_nlcd_data_centroid <- ak_nlcd_data_centroid[which(unlist(st_intersects(ak_nlcd_data_centroid, ak_hi_boundary)) == 1),] 

plot(st_geometry(ak_hi_boundary[ak_hi_boundary$STUSPS == 'AK',]))
plot(ak_nlcd_data_centroid[1:1000, 'GRIDCODE'], add = TRUE)
#ak_nlcd_data_sample <- st_intersection(ak_nlcd_data_sample, ak_hi_boundary)

# only select grid cell in AK or HI, takes hours to run
ak_nlcd_data_centroid <- ak_nlcd_data_centroid %>% select(ID, GRIDCODE, STATEFP, area, geometry)

# assign NLCD data to ak, hi tracts
ak_hi_tracts <- us_census_tract %>% dplyr::filter(STATEFP %in% c('02', '15'))
plot(st_geometry(ak_hi_tracts))

ak_nlcd_data_centroid <- st_join(ak_nlcd_data_centroid, ak_hi_tracts,
                              join = st_nearest_feature, left = TRUE)

ak_nlcd_data_centroid_df <- ak_nlcd_data_centroid %>% st_drop_geometry()
write.csv(ak_nlcd_data_centroid_df, file.path('Land_use', cleandir, 'assigned_nlcd_to_tract_ak_hi.csv'))

nlcd_data_centroid_sample <- head(ak_nlcd_data_centroid_df, 1000)          
selected_tracts <- unique(nlcd_data_centroid_sample$GEOID)
selected_tracts_sf <- ak_hi_tracts[ak_hi_tracts$GEOID %in% selected_tracts,]
plot(st_geometry(selected_tracts_sf))
plot(nlcd_data_centroid_sample[, 'GEOID'], add= TRUE)

# nlcd_data_by_tract <- nlcd_data_centroid_df %>% 
#   dplyr::group_by(GEOID, GRIDCODE) %>% 
#   dplyr::summarise(area = sum(area)) %>% 
#   tidyr::spread(GRIDCODE, area)

nlcd_data_centroid_sample_df <- nlcd_data_centroid_sample %>% st_drop_geometry()

ak_nlcd_data_centroid_df <- ak_nlcd_data_centroid_df %>% 
  dplyr::mutate(area = as.numeric(area)) # unit in m^2

ak_nlcd_data_by_tract<- ak_nlcd_data_centroid_df %>%
  tidyr::pivot_wider(names_from = GRIDCODE, values_from = area, 
                     values_fn = sum, values_fill = 0)

write.csv(ak_nlcd_data_by_tract, file.path('Land_use', cleandir, 'tract_level_land_use_ak_hi.csv'))

# check final results
us_census_tract_no_nlcd <- us_census_tract %>% 
  dplyr::filter(!GEOID %in%  unique(nlcd_data_centroid_df$GEOID))
plot(st_geometry(us_census_tract_no_nlcd))
