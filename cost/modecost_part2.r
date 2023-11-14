# Created by Michelle Chen as part of the "modecost.do" that generates "modecost.csv" on 3/4/2023
# Modified directory path by Hung-Chia Yang on 3/7/2023


mywd <- "C:/FHWA/For FHWA folks/user_and_externality_cost"
setwd(mywd)
library(dplyr)

#load("./CleanData/monetary_cost.RData")
load("./CleanData/travelmode.test.data.national_transgeo.RData")

geomicrotype <- TravelMode %>% distinct(trip_indx, o_microtype, o_geotype)
monetary_trip <- read.csv(file.path('./CleanData/trips_veh.csv')) %>%
  merge(geomicrotype, by="trip_indx")
modecost <- monetary_trip %>% group_by(o_microtype, o_geotype) %>%
  summarise(pricepermile=mean(2.41/mpg))
modecost <- modecost %>% filter(is.na(o_microtype)==F, is.na(o_geotype)==F)
modecost <- modecost %>% mutate(Microtype=paste0(o_geotype,"_", o_microtype))

write.csv(modecost, file.path("./CleanData/modecost_auto_pricepermile.csv"), row.names = F)