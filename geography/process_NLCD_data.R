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
library(raster) # for HI ige file
library(stars)

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
ak_nlcd_data_sample <- head(ak_nlcd_data_centroid, 1000)
plot(ak_nlcd_data_sample[, 'GRIDCODE'])
# get ak boundary
states <- tigris::states(cb = TRUE)
ak_boundary <- states %>% dplyr::filter(STUSPS == 'AK')
plot(st_geometry(ak_boundary))
ak_boundary <- st_transform(ak_boundary, 4326) # re-project data

# using 'st_intersects', which gives a boolean for if points within polygon
ak_nlcd_data_centroid <- ak_nlcd_data_centroid[which(unlist(st_intersects(ak_nlcd_data_centroid, ak_boundary)) == 1),] 

plot(st_geometry(ak_boundary[ak_boundary$STUSPS == 'AK',]))
plot(ak_nlcd_data_centroid[1:1000, 'GRIDCODE'], add = TRUE)
#ak_nlcd_data_sample <- st_intersection(ak_nlcd_data_sample, ak_hi_boundary)

# only select grid cell in AK or HI, takes hours to run
ak_nlcd_data_centroid <- ak_nlcd_data_centroid %>% select(ID, GRIDCODE, STATEFP, area, geometry)

# assign NLCD data to ak, hi tracts
ak_tracts <- us_census_tract %>% dplyr::filter(STATEFP %in% c('02'))
plot(st_geometry(ak_tracts))

ak_nlcd_data_centroid <- st_join(ak_nlcd_data_centroid, ak_tracts,
                              join = st_nearest_feature, left = TRUE)

ak_nlcd_data_centroid_df <- ak_nlcd_data_centroid %>% st_drop_geometry()
write.csv(ak_nlcd_data_centroid_df, file.path('Land_use', cleandir, 'assigned_nlcd_to_tract_ak.csv'))

nlcd_data_centroid_sample <- head(ak_nlcd_data_centroid, 1000)          
selected_tracts <- unique(nlcd_data_centroid_sample$GEOID)
selected_tracts_sf <- ak_tracts[ak_tracts$GEOID %in% selected_tracts,]
plot(st_geometry(selected_tracts_sf))
plot(nlcd_data_centroid_sample[, 'GEOID'], add= TRUE)

# nlcd_data_by_tract <- nlcd_data_centroid_df %>% 
#   dplyr::group_by(GEOID, GRIDCODE) %>% 
#   dplyr::summarise(area = sum(area)) %>% 
#   tidyr::spread(GRIDCODE, area)

nlcd_data_centroid_sample_df <- nlcd_data_centroid_sample %>% st_drop_geometry()

ak_nlcd_data_centroid_df <- ak_nlcd_data_centroid_df %>% 
  dplyr::mutate(area = as.numeric(area)) %>% # unit in m^2
 dplyr::select(-ID) # unit in m^2

ak_nlcd_data_by_tract<- ak_nlcd_data_centroid_df %>%
  tidyr::pivot_wider(names_from = GRIDCODE, values_from = area, 
                     values_fn = sum, values_fill = 0)

write.csv(ak_nlcd_data_by_tract, file.path('Land_use', cleandir, 'tract_level_land_use_ak.csv'))


# process Hawaii land cover
hi_nlcd_path <- 'hi_hawaii_2010_ccap_hr_land_cover20150120/hi_hawaii_2010_ccap_hr_land_cover20150120.img'
# hi_nlcd_data <- raster(file.path('Land_use', datadir, hi_nlcd_path))
# plot(hi_nlcd_data)

#tif = system.file(file.path('Land_use', datadir, hi_nlcd_path), package = "stars")
hi_nlcd_data = read_stars(file.path('Land_use', datadir, hi_nlcd_path))
hi_nlcd_data_sf <- st_as_sf(hi_nlcd_data, as_points = FALSE, merge = TRUE)
st_write(hi_nlcd_data_sf, file.path('Land_use', datadir, 'HI_NLCD_2010_P1.geojson'))
hi_nlcd_data_sample <- head(hi_nlcd_data_sf, 1000)
plot(hi_nlcd_data_sample)

hi_nlcd_path <- 'hi_hawaii_2010_ccap_hr_land_cover20150120/hi_hawaii_2010_ccap_hr_land_cover20150120.img'
# hi_nlcd_data <- raster(file.path('Land_use', datadir, hi_nlcd_path))
# plot(hi_nlcd_data)

#tif = system.file(file.path('Land_use', datadir, hi_nlcd_path), package = "stars")
hi_nlcd_path_p2 <- 'hi_hawaii_2010_ccap_hr_land_cover20150120/hi_oahu_2011_ccap_hr_land_cover20140619.img'
hi_nlcd_data_p2 = read_stars(file.path('Land_use', datadir, hi_nlcd_path_p2))
hi_nlcd_data_p2_sf <- st_as_sf(hi_nlcd_data_p2, as_points = FALSE, merge = TRUE)
hi_nlcd_data_sample <- head(hi_nlcd_data_p2_sf, 1000)
plot(hi_nlcd_data_sample)
st_write(hi_nlcd_data_p2_sf, file.path('Land_use', datadir, 'HI_NLCD_2010_oahu.geojson'))

hi_nlcd_path_p3 <- 'hi_hawaii_2010_ccap_hr_land_cover20150120/hi_niihau_2010_ccap_hr_land_cover20140930.img'
hi_nlcd_data_p3 = read_stars(file.path('Land_use', datadir, hi_nlcd_path_p3))
hi_nlcd_data_p3_sf <- st_as_sf(hi_nlcd_data_p3, as_points = FALSE, merge = TRUE)
hi_nlcd_data_sample <- head(hi_nlcd_data_p3_sf, 1000)
plot(hi_nlcd_data_sample)
st_write(hi_nlcd_data_p3_sf, file.path('Land_use', datadir, 'HI_NLCD_2010_niihau.geojson'))

hi_nlcd_path_p4 <- 'hi_hawaii_2010_ccap_hr_land_cover20150120/hi_molokai_2010_ccap_hr_land_cover20150102.img'
hi_nlcd_data_p4 = read_stars(file.path('Land_use', datadir, hi_nlcd_path_p4))
hi_nlcd_data_p4_sf <- st_as_sf(hi_nlcd_data_p4, as_points = FALSE, merge = TRUE)
hi_nlcd_data_sample <- head(hi_nlcd_data_p4_sf, 1000)
plot(hi_nlcd_data_sample)
st_write(hi_nlcd_data_p4_sf, file.path('Land_use', datadir, 'HI_NLCD_2010_molokai.geojson'))

hi_nlcd_path_p5 <- 'hi_hawaii_2010_ccap_hr_land_cover20150120/hi_maui_2010_ccap_hr_land_cover_20150213.img'
hi_nlcd_data_p5 = read_stars(file.path('Land_use', datadir, hi_nlcd_path_p5))
hi_nlcd_data_p5_sf <- st_as_sf(hi_nlcd_data_p5, as_points = FALSE, merge = TRUE)
hi_nlcd_data_sample <- head(hi_nlcd_data_p5_sf, 1000)
plot(hi_nlcd_data_sample)
st_write(hi_nlcd_data_p5_sf, file.path('Land_use', datadir, 'HI_NLCD_2010_maui.geojson'))

hi_nlcd_path_p6 <- 'hi_hawaii_2010_ccap_hr_land_cover20150120/hi_lanai_2011_ccap_hr_land_cover_20141204.img'
hi_nlcd_data_p6 = read_stars(file.path('Land_use', datadir, hi_nlcd_path_p6))
hi_nlcd_data_p6_sf <- st_as_sf(hi_nlcd_data_p6, as_points = FALSE, merge = TRUE)
hi_nlcd_data_sample <- head(hi_nlcd_data_p6_sf, 1000)
plot(hi_nlcd_data_sample)
st_write(hi_nlcd_data_p6_sf, file.path('Land_use', datadir, 'HI_NLCD_2010_lanai.geojson'))

hi_nlcd_path_p7 <- 'hi_hawaii_2010_ccap_hr_land_cover20150120/hi_kauai_2010_ccap_hr_land_cover20140929.img'
hi_nlcd_data_p7 = read_stars(file.path('Land_use', datadir, hi_nlcd_path_p7))
hi_nlcd_data_p7_sf <- st_as_sf(hi_nlcd_data_p7, as_points = FALSE, merge = TRUE)
hi_nlcd_data_sample <- head(hi_nlcd_data_p7_sf, 1000)
plot(hi_nlcd_data_sample)
st_write(hi_nlcd_data_p7_sf, file.path('Land_use', datadir, 'HI_NLCD_2010_kauai.geojson'))
# check final results
# us_census_tract_no_nlcd <- us_census_tract %>% 
#   dplyr::filter(!GEOID %in%  unique(nlcd_data_centroid_df$GEOID))
# plot(st_geometry(us_census_tract_no_nlcd))
