# Census tract shapefile compile
# updated on 11/16/23 by Xiaodan Xu, create micro-geotype boundary based on 2020 census boundary
######################################################################
# LJ add: directories so that the files are saved in the right places
mywd <- "C:/FHWA_R2/spatial_boundary"
setwd(mywd)

datadir <- "RawData"
cleandir <- "CleanData"

library(rgdal)
require(RCurl)
library(purrr) # LJ add: to use reduce
require(parallel)
# LJ add: if taRifx.geo not installed, install it once
# if ( "taRifx.geo" %in% rownames(installed.packages()) == FALSE) {
#  install.packages("remotes")
#   remotes::install_github("gsk3/taRifx.geo")
# }


#require(taRifx.geo) 
#Xiaodan notes, this package is no longer available from CRAN
# IF you need to install it, download the following GZ files and install them manually
# https://cran.r-project.org/src/contrib/Archive/taRifx.geo/
# https://cran.r-project.org/src/contrib/Archive/taRifx/

library(tidycensus) # LJ add: to get tract geometry
library(tidyverse)
library(sf)
library(sp)
library(rgeos)
library(maptools)
library(maps)
library(tigris)
library(broom) # to use tidy
library(stringr)

#setwd("~/Box/FHWA/Data/RawData/")

##########################
# CENSUS TRACT COMPILE
#######################################
#NOTE: TIGRIS files produce on 72,837 tracts for both 2017 and 2018 
# tidycensus produces 73,056
# importing and combining US census tracts
us <- unique(fips_codes$state)[1:51]
year_selected = 2020
#tidycensus
census_api_key("e74b4d8c97989e07245040ac84168a638247af9a", install = TRUE, overwrite = TRUE)
readRenviron("~/.Renviron")
options(tigris_use_cache = TRUE)

tracts <- reduce(
  purrr::map(us, function(x) {
    tracts(state = x, year = year_selected) # Xiaodan's note: using tigris package for this:https://rdrr.io/cran/tigris/man/tracts.html
    # get_acs(geography = "tract", variables = "B01003_001", 
    #         state = x, geometry = T, year = year_selected)
  }), 
  rbind
)
# tracts <- tracts %>% 
#   select(-variable, -estimate, -moe)
# exporting combined shape file. LJ add path to save the data
file_name = paste0('combined_tracts_', year_selected, '.geojson')
st_write(obj=tracts, dsn=file.path(cleandir, file_name))

tracts_df <- tracts %>% st_drop_geometry()
tracts_df <- tracts_df %>% mutate(pct_water = AWATER/(AWATER+ALAND))
class(tracts_df)
file_name = paste0('combined_tracts_', year_selected, '.csv')
write.csv(tracts_df, file.path(cleandir, file_name))
################
# COUNTY BOUNDARIES
#########################
# Natalie's original code
# combined_counties <- rbind_tigris(
#   lapply(states, function(x) {
#     counties(x, cb = TRUE, year = 2018)
#   })
# )

# LJ add:
combined_counties <- rbind_tigris(
  lapply(us, function(x) {
    counties(x, cb = TRUE, year = year_selected)
  })
)

file_name_ct = paste0('combined_county_', year_selected, '.geojson')
st_write(obj=combined_counties, dsn=file.path(cleandir, file_name_ct))

#XXu add
cbsa <- core_based_statistical_areas(cb=FALSE,resolution="500k",year=year_selected)
file_name_cbsa = paste0('combined_cbsa_', year_selected, '.geojson')
st_write(obj=file_name_cbsa, dsn=file.path(cleandir, file_name_cbsa))
# census tract 2020 - 2010 crosswalk
# link: https://www2.census.gov/geo/docs/maps-data/data/rel2020/tract/tab20_tract20_tract10_natl.txt
census_tract_crosswalk <- read.table(file.path(datadir, 
                                               'www2.census.gov_geo_docs_maps-data_data_rel2020_tract_tab20_tract20_tract10_natl.txt'), 
                                     header = TRUE, sep = "|")
write.csv(census_tract_crosswalk, file.path(cleandir, 'census_tract_crosswalk_2010_2020.csv'))

# generate GEOTYPE id and shapefile
tract_county_cbsa_crosswalk <- read.csv(file.path(cleandir, 'cleaned_lodes8_crosswalk.csv'))
tract_county_cbsa_crosswalk <- tract_county_cbsa_crosswalk %>% 
  mutate(spatial_id = ifelse(cbsa == 99999, cty, cbsa))

write.csv(tract_county_cbsa_crosswalk, file.path(cleandir, 'cleaned_lodes8_crosswalk_with_ID.csv'))

# process county
crosswalk_county_part <- tract_county_cbsa_crosswalk %>%
  filter(cbsa == 99999)
crosswalk_county_part <- crosswalk_county_part %>% mutate(is_cbsa = 0)
crosswalk_county_part <- crosswalk_county_part %>% select(spatial_id, is_cbsa) %>% distinct()
crosswalk_county_part <- crosswalk_county_part %>% mutate(spatial_id = str_pad(spatial_id, 5, pad = "0"))
# append county geometry
crosswalk_county_part <- merge(crosswalk_county_part, combined_counties,
                               by.x = 'spatial_id', by.y = 'GEOID', all.x = TRUE, all.y=FALSE) 
crosswalk_county_part <- crosswalk_county_part %>% select(spatial_id, is_cbsa, ALAND, AWATER, geometry)  

# process cbsa 
crosswalk_cbsa_part <- tract_county_cbsa_crosswalk %>% 
  filter(cbsa != 99999)
crosswalk_cbsa_part <- crosswalk_cbsa_part %>% mutate(is_cbsa = 1)
crosswalk_cbsa_part <- crosswalk_cbsa_part %>% select(spatial_id, is_cbsa) %>% distinct()
crosswalk_cbsa_part <- crosswalk_cbsa_part %>% mutate(spatial_id = str_pad(spatial_id, 5, pad = "0"))
# append county geometry
crosswalk_cbsa_part <- merge(crosswalk_cbsa_part, cbsa,
                             by.x = 'spatial_id', by.y = 'GEOID', all.x = TRUE, all.y=FALSE) 
crosswalk_cbsa_part <- crosswalk_cbsa_part %>% select(spatial_id, is_cbsa, ALAND, AWATER, geometry)  

# combine county and cbsa
combined_geotype_spatial_map <- rbind(crosswalk_county_part, crosswalk_cbsa_part)
combined_geotype_spatial_map <- st_as_sf(combined_geotype_spatial_map)

file_name_geotype = paste0('combined_geotype_unit_', year_selected, '.geojson')
st_write(obj=combined_geotype_spatial_map, dsn=file.path(cleandir, file_name_geotype))

combined_geotype_spatial_df <- combined_geotype_spatial_map %>% st_drop_geometry()
combined_geotype_spatial_df <- combined_geotype_spatial_df %>% mutate(pct_water = AWATER/(AWATER+ALAND))
class(combined_geotype_spatial_df)
file_name = paste0('combined_geotype_unit_', year_selected, '.csv')
write.csv(combined_geotype_spatial_df, file.path(cleandir, file_name))
# exporting combined shape file. LJ add path to save the data
#writeOGR(obj=combined_counties, dsn=file.path(cleandir,"Counties"), layer="combined_counties", driver="ESRI Shapefile")

# Xiaodan's notes: we didn't use the following output, please skip these lines
# combined_counties@data$GEOID <- as.character(combined_counties@data$GEOID) # convert GEOID to character
# 
# ggcbg3<-tidy(combined_counties, region = "GEOID")  # convert polygons to data.frame
# 
# fwrite(ggcbg3, file = file.path(cleandir,"Counties","all_counties_2018.csv"), row.names = F)





