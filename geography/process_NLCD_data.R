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
ak_nlcd_path <- 'emiss_shp2017/NLCD/CONUS_AK_NLCD_2011_500m_WGS.shp'
ak_nlcd_data <- st_read(file.path('Land_use', datadir, ak_nlcd_path))

ak_nlcd_data$area <- st_area(ak_nlcd_data)
ak_nlcd_data_centroid <- st_centroid(ak_nlcd_data)
ak_nlcd_data_sample <- head(ak_nlcd_data, 1000)
plot(ak_nlcd_data_sample[, 'GRIDCODE'])
# get ak boundary
states <- tigris::states(cb = TRUE)
ak_boundary <- states %>% dplyr::filter(STUSPS == 'AK')
plot(st_geometry(ak_boundary))
ak_boundary <- st_transform(ak_boundary, 4326) # re-project data
ak_nlcd_data_centroid <- st_intersection(ak_nlcd_data_centroid, ak_boundary) # only select grid cell in AK


# check final results
us_census_tract_no_nlcd <- us_census_tract %>% 
  dplyr::filter(!GEOID %in%  unique(nlcd_data_centroid_df$GEOID))
plot(st_geometry(us_census_tract_no_nlcd))
