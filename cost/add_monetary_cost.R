# Add monetary cost
# public transit fare
# Updated "trips_veh" portion of the codes by Hung-Chia Yang on 2023/06/29
rm(list = ls())
# setting working directory
mywd <- "C:/FHWA/For FHWA folks/user_and_externality_cost/"
setwd(mywd)

#figuredir <- "./Figures"
datadir <- "./RawData"
cleandir <- "./CleanData"
#tabledir <- "./Tables"


# load packages
library('pacman')
p_load(utils, foreign, reshape2, mlogit, graphics, plyr, ordinal, mltest, haven, stargazer, tidyverse, readxl)

########
#LOAD DATA
##########
trip_tract <- read.csv(file.path(datadir, './NHTS-Tract/tripct.csv'))
trips <- read.csv(file.path(datadir, './NHTS-Public/trippub.csv')) 
fare_tract <- read_excel(file.path(cleandir, 'transitfare_tract_final.xlsx'))
#fare_tract <- read_excel(file.path(cleandir, 'transitfare_tract_final.xlsx', col_types = c("text","text","numeric","numeric","numeric")))
vehicles <- read.csv(file.path(datadir, './NHTS-Public/vehpub.csv')) %>% 
  unite('veh_indx', HOUSEID:VEHID, sep = '_', remove = F)
veh_mpg <- read.csv(file.path(cleandir, 'veh_mpg.csv')) %>% 
  select(makecode, modelcode, vehtype, mpg)

#Household info from the NHTS
households <- read.csv(file.path(datadir, './NHTS-Public/hhpub.csv'))
house_tract <- read.csv(file.path(datadir,'./NHTS-Tract/hhct.csv'))

#NOTE JUNE 8 2022: these are the OLD geotypes! need to check with michelle
micro_geo <- read.csv(file.path(datadir, 'tract_cluster_results_labeled.csv')) %>%
  mutate(tract = as.character(GEOID))

#Uber fares
uber <- read_excel(file.path(datadir, 'uberfare.xlsx'))

#####################load the trip tract transit fare data ###############################
# fare_tract <- read_excel('Michelle/Data/Raw/transitfare_tract.xlsx')
# fare_microgeo <- read_excel('Michelle/Data/Raw/transitfare_micro_geo.xls')
# fare_micro <- read_excel('Michelle/Data/Raw/transitfare_micro.xls')
# 
# ###replace the NA in fare tract
# names(fare_tract)[8] <- "fare1"
# names(fare_microgeo)[4] <- "fare2"
# names(fare_micro)[3] <- "fare3"
# fare_tract <- fare_tract %>% merge(fare_microgeo, by=c("microtype", "geotype", "mode")) %>% merge(fare_micro, by=c("microtype", "mode"))
# fare_tract$fare <- fare_tract$fare1
# fare_tract$fare <- ifelse(is.na(fare_tract$fare), fare_tract$fare2, fare_tract$fare)
# fare_tract$fare <- ifelse(is.na(fare_tract$fare), fare_tract$fare3, fare_tract$fare)
# fare_tract[c("fare1", "fare2", "fare3", "microtype", "geotype", "spatial_id", "mcr_bid", "cbsa")] <- list(NULL)
# 
# rm(fare_micro, fare_microgeo, fare_tract)

####################add monetary cost for vehicle trips#########################
summary(trips$GASPRICE)

# prepare for fuel cost
trips <- trips %>% 
  unite('veh_indx', c("HOUSEID", "VEHID"), sep = '_', remove = F) %>%
  unite('trip_indx', c("HOUSEID", "PERSONID", "TDTRPNUM"), sep = '_', remove = F)

veh <- trips["veh_indx"] %>% 
  distinct() %>% 
  merge(vehicles, by="veh_indx") %>% 
  dplyr::select(MAKE, MODEL, VEHTYPE)

# prepare for parking cost
trip_tract$d_tract <- do.call(paste0, trip_tract[c("DEST_ST", "DEST_CNTY", "DEST_CT")])
trip_tract$o_tract <- do.call(paste0, trip_tract[c("ORIG_ST", "ORIG_CNTY", "ORIG_CT")])
trip_tract <- trip_tract %>% unite('trip_indx', HOUSEID:TDTRPNUM, sep = '_', remove = F)
trips <- trips %>% merge(trip_tract, by="trip_indx")
parking <- trips %>% dplyr::select(trip_indx, d_tract, WHYTRP1S)
# write.csv(parking, "Michelle/Data/Raw/parking.csv")

# process mpg data
trips <- trips %>% 
  left_join(vehicles, by="veh_indx")
trips_veh <- trips %>% 
  dplyr::select(trip_indx, veh_indx, MAKE, MODEL, VEHTYPE) %>% 
  merge(trip_tract, by="trip_indx")

names(veh_mpg) <- c("MAKE", "MODEL", "VEHTYPE", "mpg")
trips_veh <- trips_veh %>% 
  left_join(veh_mpg, by=c("MAKE", "MODEL", "VEHTYPE"))

house_tract <- house_tract %>% mutate(tract=HHSTFIPS*10^9+HHCNTYFP*10^6+HHCT)
house_tract$tract <- as.character(house_tract$tract)
house_tract$tract[house_tract$HHSTFIPS<10] <- paste("0", house_tract$tract[house_tract$HHSTFIPS<10], sep="")
house_tract <- house_tract %>% dplyr::select(HOUSEID, tract)
trips_veh <- trips_veh %>% merge(house_tract, by="HOUSEID")
micro_geo$tract[micro_geo$fips_st<10] <- paste("0", micro_geo$tract[micro_geo$fips_st<10], sep="")
micro_geo <- micro_geo %>% dplyr::select(tract, geo_micro_type, microtype, geotype)
trips_veh <- trips_veh %>% merge(micro_geo, by="tract")
mpg_microgeo <- trips_veh %>% 
  dplyr::group_by(geo_micro_type) %>%
  dplyr::summarise(mpg1=mean(mpg, na.rm=T))
trips_veh <- trips_veh %>% merge(mpg_microgeo, by="geo_micro_type")
trips_veh$mpg[is.na(trips_veh$mpg)] <- trips_veh$mpg1[is.na(trips_veh$mpg)]
trips_veh <- trips_veh %>% dplyr::select(trip_indx, mpg)

write.csv(trips_veh, file.path(cleandir,"trips_veh.csv"), row.names = F)

# load cleaned parking data
# parking <- na.omit(parking)
# parking$d_tract <- as.numeric(parking$d_tract)
# parking$state <- floor(parking$d_tract/(10^9))
# parking$d_tract <- as.character(parking$d_tract)
# parking$d_tract[parking$state<10] <- paste("0", parking$d_tract[parking$state<10], sep="")
# trips_parking <- trips %>% dplyr::select(trip_indx, d_tract, WHYTRP1S)
# names(trips_parking)[3] <- "whytrp1s"
# names(parking)[3] <- "whytrp1s"
# trips_parking <- trips_parking %>% merge(parking, 
#                                          by=c("d_tract", "whytrp1s"), all.x=TRUE)
# trips_parking$parking[is.na(trips_parking$parking)] <- 0
# trips_parking <- trips_parking %>% dplyr::select(trip_indx, parking)
# 
# # load cleaned public transit fare data
# trips_fare <- trips %>% 
#   dplyr::select(trip_indx, o_tract) %>% 
#   merge(fare_tract, by.x="o_tract", by.y="tractcode")
# 
# # load uber fare
# trips_uber <- trips %>% dplyr::select(trip_indx, o_tract)
# uber <- uber %>% dplyr::select(tractcode, minimumfare, priceperminute, pricepermile)
# trips_uber <- trips_uber %>% 
#   merge(uber, by.x="o_tract", by.y="tractcode") %>%
#   mutate(mode = "taxi")
# 
# trips_veh <- trips_veh %>% 
#   merge(trips_parking, by="trip_indx") %>%
#   mutate(mode = "hv")
# 
# trips_bike <- trips %>% 
#   dplyr::select(trip_indx) %>%
#   mutate(bikecost = 0.585,
#          mode = "bike")
# 
# trips_walk <- trips %>% 
#   dplyr::select(trip_indx) %>%
#   mutate(walkcost = 0,
#          mode = "walk")
# 
# save(trips_fare, trips_veh, trips_uber, trips_bike, trips_walk, file = file.path(cleandir, 'monetary_cost.RData'))

#load("Michelle/Data/Clean/monetary_cost.RData")
# reformat the data
#trip_indx <- trips_uber["trip_indx"]
#trips_veh <- trips_veh %>% merge(trip_indx, by="trip_indx")
#trips_bike <- trips_bike %>% merge(trip_indx, by="trip_indx")
#trips_walk <- trips_walk %>% merge(trip_indx, by="trip_indx")

#monetary <- bind_rows(trips_fare, trips_uber, trips_veh, trips_bike, trips_walk)
#save(monetary, file = file.path(outdir, 'monetary_trips.RData'))
#load(file = file.path(outdir, 'monetary_trips.RData'))
#rm(list=setdiff(ls(), c("monetary", "datdir", "outdir", "tabledir")))

########################add monetary cost to existing mode choice data############################
# load("C:/FHWA/For FHWA folks/CleanData/TravelMode.modmat.National.small_transgeo.RData")
# TravelMode <- TravelMode.modmat
# rm(TravelMode.modmat)
# trips_fuel <- trips %>% dplyr::select(trip_indx, GASPRICE)
# 
# TravelMode <- TravelMode %>% 
#   merge(monetary, by=c("trip_indx", "mode")) %>% 
#   merge(trips_fuel, by=c("trip_indx")) %>% 
#   mutate(veh_cost=parking+(trpmiles/mpg)*GASPRICE/100,
#          taxi_cost=inv_time*priceperminute+trpmiles*pricepermile)
# 
# TravelMode$taxi_cost[TravelMode$taxi_cost<TravelMode$minimumfare & TravelMode$mode=="taxi"]<- TravelMode$minimumfare[TravelMode$taxi_cost<TravelMode$minimumfare  & TravelMode$mode=="taxi"]
# TravelMode <- TravelMode %>% 
#   mutate(bike_cost=trpmiles*0.585, 
#          transit_cost=fare)
# TravelMode$transit_cost[TravelMode$trpmiles>=10 & (TravelMode$mode=="rail_l" | TravelMode$mode=="rail_c")]<- 
#   TravelMode$maxi[TravelMode$trpmiles>=10 & (TravelMode$mode=="rail_l" | TravelMode$mode=="rail_c")]
# 
# TravelMode <- TravelMode %>% mutate(cost=case_when(mode=="hv"~veh_cost,
#                                                    mode=="taxi"~taxi_cost, 
#                                                    mode=="bike"~bike_cost, 
#                                                    mode=="bus"~fare, 
#                                                    mode=="rail_l"~transit_cost, 
#                                                    mode=="rail_c"~transit_cost,
#                                                    TRUE~0))
# TravelMode$cost[is.na(TravelMode$cost)]<-TravelMode$fare[is.na(TravelMode$cost)]
# summary(TravelMode$cost)
# TravelMode <- TravelMode %>% dplyr::select(-fare:-transit_cost)
# 
# ## Add number of household members on the trip
# trips <- trips %>% dplyr::select(trip_indx, TRPHHACC)
# TravelMode <- TravelMode %>% merge(trips, by="trip_indx")
# TravelMode$TRPHHACC[TravelMode$TRPHHACC==-9] <- 0
# TravelMode$TRPHHACC <- TravelMode$TRPHHACC+1
# TravelMode <- TravelMode %>% mutate(inv_time2=inv_time*TRPHHACC, access_time2=access_time*TRPHHACC, wait_time2=wait_time*TRPHHACC, cost2=cost)
# TravelMode$cost2[TravelMode$mode=="rail_c"] <- TravelMode$cost[TravelMode$mode=="rail_c"]*TravelMode$TRPHHACC[TravelMode$mode=="rail_c"]
# TravelMode$cost2[TravelMode$mode=="rail_l"] <- TravelMode$cost[TravelMode$mode=="rail_l"]*TravelMode$TRPHHACC[TravelMode$mode=="rail_l"]
# TravelMode$cost2[TravelMode$mode=="bus"] <- TravelMode$cost[TravelMode$mode=="bus"]*TravelMode$TRPHHACC[TravelMode$mode=="bus"]
# 
# table(TravelMode$TRPHHACC)
# TravelMode <- TravelMode %>% mutate(person1 = case_when(TRPHHACC==1~1, TRUE~0),
#                                     person2more = case_when(TRPHHACC>=2~1, TRUE~0),
#                                     person24 = case_when(TRPHHACC>=2 & TRPHHACC<=4 ~1, TRUE~0),
#                                     person5more = case_when(TRPHHACC>=5~1, TRUE~0))
# ## change the cost bike into 0
# TravelMode$cost[TravelMode$mode=="bike"]<-0
# TravelMode$cost2[TravelMode$mode=="bike"]<-0
# save(TravelMode, file = file.path(cleandir, 'TravelMode_full.RData'))
# 
# ##### graphs for distribution
# p1 <- ggplot(trips_veh, aes(x=mpg)) + 
#   geom_histogram(aes(y=..density..), binwidth=2, color="darkblue", fill="lightblue") +
#   labs(title="Distribution of vehicles' MPG for all trips",x="MPG(miles/gallon)", y = "Density") +
#   theme(plot.title = element_text(hjust = 0.5))
# p1
# 
# TravelMode <- TravelMode %>% filter(mode=="taxi") %>% select(trip_indx, cost)
# trips_uber <- trips_uber %>% 
#   merge(TravelMode, by="trip_indx") %>% filter(cost<200)
# 
# p2 <- ggplot(trips_uber, aes(x=cost)) + 
#   geom_histogram(aes(y=..density..), binwidth=10, color="darkgreen", fill="lightgreen") +
#   labs(title="Distribution of taxi fares",x="Cost (dollar)", y = "Density") +
#   theme(plot.title = element_text(hjust = 0.5))
# p2
# 
# uber <- read.csv(file.path(datadir, 'Uber/Uber_fare.csv')) # Note: this is different than the other urber data- ask Michelle
# p3 <- ggplot(uber, aes(x=Price.Per.Mile)) + 
#   geom_histogram(aes(y=..density..), binwidth=0.1, color="darkgreen", fill="lightgreen") +
#   labs(x="Price per mile (dollar)", y = "Density") +
#   theme(plot.title = element_text(hjust = 0.5))
# p3
# 
# p4 <- ggplot(trips_fare %>% filter(mode=="bus"), aes(x=fare)) + 
#   geom_histogram(aes(y=..density..), binwidth=0.5, color="darkgreen", fill="lightgreen") +
#   labs(title="Distribution of bus fares",x="Cost (dollar)", y = "Density") +
#   theme(plot.title = element_text(hjust = 0.5))
# p4
