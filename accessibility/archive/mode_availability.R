rm(list=ls())
#setwd("D:/Box Sync/FHWA/Task2/")

mywd <- "C:/FHWA/For FHWA folks/opportunity_and_mode_availability"
setwd(mywd)

#figuredir <- "./Figures"
datadir <- "RawData"
cleandir <- "CleanData"

# load packages
library('pacman')
p_load(utils, foreign, pastecs, mlogit, graphics, VGAM, aod, plotrix, Zelig, Zelig, 
       vctrs, maxLik, plyr, MASS, ordinal, mltest, haven, stargazer, stringr, tidyverse)
p_load(gWidgets2, gWidgets2tcltk, miscTools, lmtest, dplyr, BiocManager)
p_load(ggplot2, scales)
### Xiaodan's notes: MASS and dplyr have functions with the same name, be careful about which one is used


######
# LOAD DATA
##########
# crosswalk data for GEOIDs --> Xiaodan's notes -- this is simply a better way to find geoid
geoid_lookup <- read.csv(file.path(cleandir, "us_xwalk_tract_2017_withID.csv"))

# Rail density
rail <- read.csv(file.path(datadir, "rail_nodes_per_tract.csv"))

# CNT data for bus service 
# Source: email from CNT (Peter Haas)
cnt <- read.csv(file.path(datadir, "us_tracts_ht_inputs.csv"))

# Distance between each bus stop --> Xiaodan's note: input not used 
# bus_shape <- read.csv(file.path(datadir, "distance_between_stops.csv"))

#bikeshare station density
bike <- read.csv(file.path(datadir, "bikedensity_area2017.csv"))

# distance to rail stations 
distance_rail <- read.csv(file.path(datadir, "distance_point.csv"))

#distance to bus stops  --> Xiaodan's note: input not used 
# distance_bus <- read.csv(file.path(cleandir, "distance_bytract.csv"))

# Subset of LA census tracts only -->Xiaodan note: not needed 
# la_tract <- read.csv(file.path(datadir, 'LA_cbsa_tracts.csv'))

#=======================================Mode Accessbility===================================#
# the visualization code are not activated
# hist(cnt$tas_trips_per_week[cnt$tas_trips_per_week<500])
# numtrips <- cnt %>% 
#   filter(tas_trips_per_week<50) %>% 
#   ggplot(aes(x=tas_trips_per_week)) + 
#   geom_histogram(binwidth = 0.5, color="black", fill="white")+
#   scale_x_continuous(limits=c(-0.000001, 50)) + 
#   scale_y_continuous(limits=c(0, 100))
# 
# numtrips <- cnt %>%
#   ggplot(aes(x=tas_trips_per_week)) + 
#   geom_histogram(binwidth = 500, color="black", fill="white") #+
#   # scale_x_continuous(limits=c(-0.000001, 50)) + scale_y_continuous(limits=c(0, 100))
# numtrips
# 
bus <- cnt %>%
  dplyr::mutate(bus=case_when(tas_trips_per_week<14 ~ 0, TRUE ~ 1)) %>%
  dplyr::select(stfid, bus)
# 
# bus_shape <- bus_shape %>% 
#   group_by(geoid) %>% 
#   summarise(nubmer=n()) %>% 
#   filter(nubmer<150) %>% 
#   ggplot(aes(x=nubmer)) + 
#   geom_histogram(binwidth = 5, color="black", fill="white")

bike <- bike %>% 
  dplyr::mutate(bike=case_when(density_land>0~1, TRUE~0)) %>% 
  dplyr::select(geoid, bike)

rail <- rail %>% 
  dplyr::group_by(geoid) %>% 
  dplyr::summarise(rail = sum(numpoints)) %>% 
  dplyr::mutate(rail=case_when(rail==0~0, TRUE~1)) %>% 
  dplyr::select(geoid, rail)

# plot <- rail %>% 
#   filter(rail<200) %>%
#   ggplot(aes(x=rail)) + 
#   geom_histogram(binwidth = 10, color="black", fill="white")
# plot

geoid_lookup <- geoid_lookup %>% 
  dplyr::select(GEOID)

test <- bike %>% 
  merge(geoid_lookup, by.x="geoid", by.y="GEOID") %>% 
  merge(rail, by.x="geoid", by.y="geoid", all.x = T)
test$rail[is.na(test$rail)] <- 0

test <- test %>% 
  merge(bus, by.x="geoid", by.y="stfid", all.x = T)
test$bus[is.na(test$bus)] <- 0

write.csv(test, file.path(cleandir,"modeaccessbility.csv"), row.names = T)

#=======================================Mode Accessbility===================================#
# rail distance
tract <- read.csv(file.path(cleandir, "contiguous_microtype_blobs_key_labeled_transgeo.csv"))
tract <- tract %>% 
  merge(rail, by.x="GEOID", by.y="geoid") %>% 
  filter(rail>0) %>% 
  merge(distance_rail, by.x="GEOID", by.y="geoid")


tract_geo <- tract %>% 
  dplyr::group_by(geotype, microtype) %>% 
  dplyr::summarise(distance=mean(distance))

write.csv(tract_geo, file.path(cleandir,"raildistance.csv"), row.names = T)


##### Xiaodan note: code below is not needed for national analysis ####
# CA rail distance
# tract_geo <- tract %>% 
#   mutate(state=floor(GEOID/10^9)) %>% 
#   filter(state==6) %>% 
#   group_by(geotype, microtype) %>% 
#   summarise(distance=mean(distance))
# 
# # bus distance
# tract <- read.csv(file.path(cleandir,"contiguous_microtype_blobs_key_labeled_transgeo.csv"))
# tract = tract %>% 
#   merge(bus, by.x="GEOID", by.y="stfid") %>% 
#   filter(bus>0) %>% 
#   merge(distance_bus, by.x="GEOID", by.y="geoid")
# 
# tract_geo <- tract %>% 
#   group_by(geotype, microtype) %>% 
#   summarise(distance=mean(distance))
# 
# # CA bus distance
# tract_geo <- tract %>% 
#   mutate(state=floor(GEOID/10^9)) %>% 
#   filter(state==6) %>% 
#   group_by(geotype, microtype) %>% 
#   summarise(distance=mean(distance))
# 
# # la bus distance
# tract_geo <- tract_la %>%
#   mutate(state=floor(GEOID/10^9)) %>% 
#   filter(state==6) %>% 
#   group_by(geotype, microtype) %>% 
#   summarise(distance=mean(distance))
# 
# # LA
# tract_la <- tract %>% 
#   merge(la_tract, by.x="GEOID", by.y="tract")
# 
# tract_geo_la <- tract_la %>% 
#   group_by(geotype, microtype) %>% 
#   summarise(distance=mean(distance))
