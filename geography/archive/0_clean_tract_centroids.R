# Compute geometric centroid since a bunch of the population centroids are not listed in the census website

mywd <- "C:/FHWA/For FHWA folks"
setwd(mywd)

#LJ add: fixed the paths
#figuredir <- "./Figures"
datadir <- "RawData"
cleandir <- "CleanData"
#tabdir <- './Tables'
#resultsdir <- './Results'

library(tidycensus)
library(data.table)
library(sf)
library(tidyverse)

########
# LOAD DATA
#############
#census tract boundaries
census_api_key("e74b4d8c97989e07245040ac84168a638247af9a", overwrite = TRUE)
readRenviron("~/.Renviron")
options(tigris_use_cache = TRUE)

# get FIPS codes
us <- unique(fips_codes$state)[1:51]

tracts <- reduce(
  map(us, function(x) {
    get_acs(geography = "tract", variables = "B01003_001", 
            state = x, geometry = TRUE, year = 2017)
  }), 
  rbind
)

# compute geoemtric centroids for those that are missing the population centroids
# change to projected CRS
tracts2 = tracts
tracts2$centroid = st_transform(tracts2, crs = 3857) %>% # projected WGS 84
  st_centroid() %>% 
  st_transform(., crs  = 4269) %>% #transform back to NAD83
  st_geometry() #assign geometry to centroids

#tracts2 = st_transform(tracts2, crs  = 4269)
tracts2 = tracts2 %>%
  select(GEOID, centroid) %>%
  mutate(LATITUDE = unlist(map(tracts2$centroid,2)),
         LONGITUDE = unlist(map(tracts2$centroid,1)),
         STATEFP = substr(GEOID,1,2), #separate GEOID to match census population files
         COUNTYFP = substr(GEOID,3,5),
        TRACTCE = substr(GEOID,6,11)) %>%
  select(-centroid, -GEOID) %>%
  st_drop_geometry() %>%
  filter(!is.na(LATITUDE)) 

write.csv(tracts2,file.path(cleandir, "tracts_centroids_geometric.csv"), row.names = F)








