#Update 10.19 combine railways: 
#Update 11.17 add disability: 
#Update 2022.1.5 CA:
#Last updated: 10/2022 by Sayeeda Ayaz
############################# Only use the trips having access to all modes#######################################
rm(list=ls())
setwd('C:/FHWA/For FHWA folks')
datadir <- './Mode_choice_estimation/Data'


# CA.TF = F


# initialization of the clustering
library('pacman')

p_load(dplyr, MASS, reshape2,cowplot,rgdal, # rgdal: reading shapefiles
       ggplot2, corrplot, raster,readxl,RColorBrewer,fmsb) # sp: reading shapefiles; RColorBrewer color paletters; fmsb: radar charts
p_load(ncdf4,rgdal,ggmap,sp,shapefiles, maps, sf, fields, Imap, raster,readxl,
       tictoc,stargazer,psych,GPArotation,spdplyr,sp,shapefiles,tmap) # tictoc: time procedures

p_load(cluster,factoextra,DandEFA, xtable,psychTools, # DandEFA: for dandelion plots of EFA factors; factoextra: metrics to choose number of clusters
       xtable,psychTools,aCRM,clusterCrit,data.table,tigris,DAAG, fastDummies,tidyr,mlogit, filling, FNN)

#########
# LOAD DATA
#########
parking <- read.csv(file=file.path(datadir, 'parking_tract.csv')) 

trips = read.csv(file.path(datadir, 'trippub.csv')) 

# Person weights from nhts
wts = read.csv(file = file.path(datadir,'perwgt.csv')) 

# train distance data merged to trip locations.
#load(file= file.path(datadir,'NHTS_tract_origin_train_dist_transgeo.RData')) 

# Taxi fare
taxi <- read.csv(file.path(datadir, 'uber_tract_matched.csv')) 

# Bikeshare station density
bikeshare <- read.csv(file.path(datadir, "bikedensity_area2017.csv")) 

# Houshold GEOIDs?
location <- read.csv(file = file.path(datadir, "hhct.csv" )) 

# Vehicles from the NHTS
vehicles <- read.csv(file=file.path(datadir, 'vehpub.csv')) 

# public transit fares
public <- read_xlsx(file.path(datadir, 'transitfare_tract_final.xlsx'), guess_max = 10000) 

trips <- trips %>% 
  unite('trip_indx', HOUSEID:TDTRPNUM, sep = '_', remove = F) # 923,572 trips


colnames(trips) = tolower(colnames(trips))
names(trips)
dim(trips)

#rename to CA convention
trips = trips %>%
  dplyr::rename(sampno = houseid, 
                perno = personid) %>% 
  filter(!is.na(sampno) & !is.na(perno) & !is.na(tdtrpnum))
dim(trips) 

# use WTPERFIN - person final weights
wts = wts %>%
  dplyr::select(HOUSEID, PERSONID,  WTPERFIN) %>% 
  dplyr::rename(sampno = HOUSEID, 
                persno = PERSONID,
                wtperfin = WTPERFIN)
##############################################
# 1. extract relevant variables
select.vars = 
  c("trip_indx",
    "sampno", # household id
    "perno",  # person id within the household
    "tdtrpnum", # trip id
    # "wttrdfin", #from tripchains dataset
    "strttime", # trip start time
    'endtime', # trip end time
    "trvlcmin", # trip duration or travel time (calculated travel time)
    "trpmiles" , # total trip length in miles, or travel distance
    "trptrans", # derived trip mode, only single mode is selected.
    "trwaittm", # transit wait time in minutes, 0-60, set other values to NA
    "tracctm" , # trip time to transit station in minutes
    "drop_prk" , # park or dropped at transit station, 1-parked, 2-dropped, 3, bike/walked, other -NA
    "tregrtm" , # time used from transit to destination in minutes [0,180]
    # "tour", #from tripchains dataset
    # "tour_seg",
    # "tourtype",
    # "trpcnt",
    # "stops",
    # "tour_flg",
    # "trip_chain", #New variable
    #    "tdwknd" , # is it a weekend trip. Note that generally there is only weekday data
    #    "vmt_mile" ,# Trip distance in miles for personally driven vehicle trips, derived from route geometry returned by Google Maps API
    "whytrp1s" # trip purpose summary
  ) 

dat = trips[,select.vars] 


#%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% merge with the trip chain data set %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
#Load the 2017 NHTS tripchain data set
# tripchains = read.csv(file.path(datadir, 'chntrp17.csv')) 
# 
# # Has the variable TRIP_CHAIN with 3 levels: 0 = 1 legged trip, 1 = 2-legged trip
# tripchainingindex = read.csv(file.path(datadir, 'Trip_chaining_index.csv')) 
# 
# tripchains = tripchains %>%
#   merge(tripchainingindex) 
# 
# # Create the tour index to perform the trip chaining
# trip_tour = tripchains %>% 
#   unite('tour_indx', HOUSEID, PERSONID, TOUR, sep = '_', remove = F)%>%
#   unite('trip_indx', HOUSEID:TDTRPNUM, sep = '_', remove = F) 

# colnames(trip_tour) = tolower(colnames(trip_tour))
# 
# trip_tour = trip_tour%>%
#   select(-houseid, -personid, -tdtrpnum)

# dat = dat %>%
#   merge(trip_tour, by = "trip_indx")

#### take a look at wait_time by different bus modes ####
waittime = dat %>%
  filter(trwaittm>0)%>%
  group_by(trptrans) %>%
  summarise(mean_wait_time = mean(trwaittm)) # mean_wait_time = 9.260688

######## only keep the city bus mode as bus mode #####
dat = dat %>%
  filter(!(trptrans %in% c(10,12,13,14))) %>% # exclude school bus, charter bus, dial a ride and intercity bus
  filter(!(trptrans == 11 & trwaittm > 50))  # exclude bus trips with more than 50 min waiting time, 10% of the bus trips,
#909,863 trips
######################################################
# 2. code up mode_chosen and upper_level_nest
dat = dat %>% 
  dplyr::mutate(mode_chosen = case_when(
    trptrans %in% c(3,4,5,6,18,9) ~ 'hv', # human driving vehicle
    trptrans %in% c(2) ~ 'bike',
    trptrans %in% c(7,8) ~ 'scooter',
    trptrans %in% c(1) ~ 'walk',
    trptrans %in% 10:14 ~ 'bus', # now we only have 11 in the data
    trptrans %in% c(15, 16) ~ 'rail',
    trptrans %in% c(17) ~ 'taxi',
    TRUE ~ 'other'
  )) %>%
  dplyr::mutate(upper_level_nest = case_when(
    mode_chosen %in% c('hv') ~ 'auto',
    mode_chosen %in% c('bike','scooter','walk') ~ 'micromobility',
    mode_chosen %in% c('bus','rail','taxi') ~ 'transit',
    TRUE ~ 'other'
  ))

dat = as.data.frame(dat)  

################################
# 3. code up trip purposes

dat = dat %>%
  dplyr::mutate(trip_purpose = case_when(
    whytrp1s %in% c(50) ~ 'social',
    whytrp1s %in% c(1) ~ 'home',
    whytrp1s %in% c(10) ~ 'work',
    whytrp1s %in% c(40) ~ 'shopping',
    whytrp1s %in% c(70) ~ 'transp_someone',
    whytrp1s %in% c(30) ~ 'medical',
    whytrp1s %in% c(80) ~ 'meals',
    whytrp1s %in% c(20) ~ 'school',
    TRUE ~ 'other'
  )) 

################################
# 4. code up travel time in minute, including total time, access time, and in-transit/vehicle time
dat = dat %>%
  filter(mode_chosen != 'other') %>% # excluding  trips in "other" mode
  filter(mode_chosen != 'scooter') %>% # not considering scooter mode for now
  filter(trvlcmin > 0) %>% # excluding  trips with no travel time obtained
  filter(upper_level_nest %in% c('auto','micromobility') | 
           (upper_level_nest == 'transit' & tracctm >= -1 & trwaittm >= -1)) %>% # drop transit trips with NAs in access and waittime.
  dplyr::mutate(total_travel_time = as.numeric(trvlcmin)) %>% # total travel time in min
  dplyr::mutate(access_time = case_when(  # travel time to transit station in minutes, if the mode is hv, bike/walk, this is 0
    #    upper_level_nest %in% c('auto','micromobility') ~ 0,                
    #    upper_level_nest == 'transit' &  tracctm == -1 ~ 0,
    upper_level_nest == 'transit' &  tracctm >0 ~ as.numeric(tracctm) ,
    TRUE ~ 0 
  )) %>%
  dplyr::mutate(wait_time = case_when( ## wait time at transit stations in minutes, if the mode is hv, bike/walk, this is 0
    upper_level_nest == 'transit' &  trwaittm >0 ~ as.numeric(trwaittm) ,
    TRUE ~ 0
  )) %>%
  dplyr::mutate(egress_time = case_when(
    upper_level_nest == 'transit' &  tregrtm >0 ~ as.numeric(tregrtm),
    TRUE ~ 0 
  ))%>%
  dplyr::mutate(inv_time = total_travel_time - wait_time - access_time,
                inv_time2 = total_travel_time - wait_time - access_time - egress_time) %>% # in transit time if travel by transit or travel time if other modes
  filter(inv_time > 0) %>%  # remove transit trips with derived in transit time less than 0     
  filter(trpmiles >0) %>% # excluding 91 trips that had NA in travel distance
  dplyr::mutate(inv_speed = as.numeric(trpmiles)/inv_time) %>%
  filter(inv_speed < 100/60) # removed 257 trips with derived speed exceedingly large (i.e. >100 mile/h),

# 896,341 trips after dropping cases with invalid travel time

###################################################
# 5. Generate travel time bins (now not decided): morning rush hour [6am, 9am), evening rush hour [4pm, 7pm], and other hours
# NOTE NDP 5/20/21: We need to preserve the variable "strttime" so we can create time bins on the fly
dat = dat %>%
  dplyr::mutate(start_time_bin = case_when(
    strttime >=600 & strttime <900 ~ 'morning_rush',
    strttime >=1600 & strttime <=1900 ~ 'evening_rush',
    TRUE ~ 'other_time'
  )) 

#######################################################################
# 6. to generate variables for alternative modes, need to derive travel speed by each mode, then summarize by micro_geotype o-d pairs and by travel start-time bin
#    then generate travel time or inv_time for each mode
master.dat = dat
colnames(master.dat) = tolower(colnames(master.dat))

load(file = file.path(datadir,'NHTS_tract_od_with_label_no_location.RData')) 

#od.tr.purpose2 = unfold.list(od.tr.purpose2)
# map to our final cluster labels (transgeo version already mapped to final labels)
#od.tr.purpose2$o_microtype = map_labels(od.tr.purpose2,raw_colname = 'o_microtype',microtype = T)
#od.tr.purpose2$d_microtype = map_labels(od.tr.purpose2,raw_colname = 'd_microtype',microtype = T)
#od.tr.purpose2$o_geotype = map_labels(od.tr.purpose2,raw_colname = 'o_geotype',microtype = F)
#od.tr.purpose2$d_geotype = map_labels(od.tr.purpose2,raw_colname = 'd_geotype',microtype = F)

tmp = od.tr.purpose2 %>% 
  unite('trip_indx', HOUSEID:TDTRPNUM, sep = '_', remove = F) %>% 
  dplyr::select(-whytrp1s) %>% 
  unite('o_geo_mic', o_geotype, o_microtype, remove = F) %>%
  unite('d_geo_mic', d_geotype, d_microtype, remove =F) %>%
  unite('od_geo_mic', o_geo_mic, d_geo_mic,remove = F) %>%
  dplyr::rename(sampno = HOUSEID, 
                perno = PERSONID,
                tdtrpnum = TDTRPNUM) %>%
  inner_join(dat) %>%   
  filter((d_microtype >0 & o_microtype >0) & !is.na(d_microtype) & !is.na(o_microtype) ) #remove NAs and water.

dat = tmp 
# prepare inv_speed (miles/min)
#dat = dat %>%
#  filter(trpmiles >0) 

# observed speed by mode, start_time, and od pairs, 
# we can see travel speed for non-auto modes are not observed in all o-d pairs of auto trips 
# so for the time being, generate average speed by mode as a function of start.time and mode
# testing = F
# if(testing){ 
#   mode.speed = dat %>%
#     group_by(start_time_bin,od_geo_mic, mode_chosen,o_microtype,d_microtype) %>%
#     summarise(mean.speed = mean(inv_speed),
#               n.sample = n())
#   
#   mode.speed %>%
#     group_by(start_time_bin, mode_chosen) %>%
#     summarise(N_unique_od_pairs = n_distinct(od_geo_mic)) %>%
#     group_by( mode_chosen) %>%
#     summarise(N_unique_od_pairs = n_distinct(od_geo_mic))
#   
#   # how about group by od-micriotype pairs --> We can see there are a number of locations with no taxi and rail speed observed
#   groupby.od.micro = F
#   if(groupby.od.micro ){
#     mode.speed = dat %>%
#       group_by(start_time_bin,o_microtype,d_microtype, mode_chosen) %>%
#       summarise(mean.speed = mean(inv_speed),
#                 n.sample = n())
#     
#     ggplot(data = mode.speed, aes(x = mode_chosen, y = mean.speed)) +
#       geom_boxplot(aes(color = start_time_bin)) +
#       theme_bw() +
#       facet_grid(o_microtype ~ d_microtype)
#   }
#   

#   # look at the mode speed (first averaged over micro-geotype pairs)
#   ggplot(data = mode.speed, aes(x = mode_chosen, y = mean.speed*60)) +
#     geom_boxplot(aes(color = start_time_bin))
#   
#   # look at the average raw speed 
#   ggplot(data = dat, aes(x = mode_chosen, y = inv_speed*60)) +
#     geom_boxplot(aes(color = start_time_bin))+
#     theme_bw() +
#     facet_grid(o_microtype ~ d_microtype)
# }



# So for the time being, generate average speed by mode 
# as a function of start.time and mode, and od microtype pairs
# if no observation at od microptype pairs, use the mean value.
# LJ 8/12/2020 change to harmonic means harmonic.mean(x)


mode.speed = dat %>%
  dplyr::group_by(start_time_bin, mode_chosen,o_microtype,d_microtype) %>%
  dplyr::summarise(mean.speed = harmonic.mean(inv_speed),
                   mean.speed.m_h = harmonic.mean(inv_speed) * 60) 

# mode.speed1 = dat %>%
#   group_by(mode_chosen,o_geotype,d_geotype) %>%
#   summarise(mean.speed = harmonic.mean(inv_speed),
#             mean.speed.m_h = harmonic.mean(inv_speed) * 60)
# 
# mode.speed2 = dat %>%
#   group_by(mode_chosen,o_geotype) %>%
#   summarise(mean.speed = harmonic.mean(inv_speed),
#             mean.speed.m_h = harmonic.mean(inv_speed) * 60)
# 
# mode.speed3 = dat %>%
#   group_by(mode_chosen,o_microtype,d_microtype) %>%
#   summarise(mean.speed = harmonic.mean(inv_speed),
#             mean.speed.m_h = harmonic.mean(inv_speed) * 60)
# 
# mode.speed4 = dat %>%
#   group_by(mode_chosen,o_microtype) %>%
#   summarise(mean.speed = harmonic.mean(inv_speed),
#             mean.speed.m_h = harmonic.mean(inv_speed) * 60)


mode.speed.mean = dat %>%
  dplyr::group_by(start_time_bin, mode_chosen) %>%
  dplyr::summarise(mean2.speed = harmonic.mean(inv_speed),
                   mean2.speed.m_h = harmonic.mean(inv_speed) * 60 ) 


# observed access.time and wait.time by origin, start.time. and mode chosen
# full set should include 6*6 (origin micro-geo types) * 3 (time bin) * 6 (modes) = 1134
# instead we observe 43 origin micro-geo types
# there are more origin types in hv than other modes, so observed wait time and access time will not cover all the origins in the auto (hv) trips
#so for the time being, use wait and access time as a function of start.time and mode. 
# testing = F
# if(testing){ 
#   transit.times = dat %>% 
#     group_by(start_time_bin,o_geo_mic, mode_chosen) %>%
#     summarise(mean.access.time = mean(access_time),
#               mean.waittime = mean(wait_time))
#   
#   transit.times %>%
#     group_by(start_time_bin,mode_chosen) %>%
#     summarise(N_origin = n_distinct(o_geo_mic)) 
#   
#   transit.times %>%
#     group_by(mode_chosen) %>%
#     summarise(N_origin = n_distinct(o_geo_mic))
#   
#   ggplot(data = transit.times, aes(x = mode_chosen, y = mean.access.time)) +
#     geom_boxplot(aes(color = start_time_bin))
#   
#   ggplot(data = transit.times, aes(x = mode_chosen, y = mean.waittime)) +
#     geom_boxplot(aes(color = start_time_bin))
#   
#   # how about only as a function of origin microtype rather than micro-geotype combinations?
#   # --> we now have everything observed. 
#   
#   groupby.o.micro = F
#   if(groupby.o.micro){
#     transit.times = dat %>% 
#       group_by(start_time_bin,o_microtype, mode_chosen) %>%
#       summarise(mean.access.time = mean(access_time),
#                 mean.waittime = mean(wait_time))
#     
#     ggplot(data = transit.times[transit.times$o_microtype!=7,], aes(x = mode_chosen, y = mean.access.time)) +
#       geom_boxplot(aes(color = start_time_bin)) +
#       theme_bw() +
#       facet_wrap(.~o_microtype)
#     
#     ggplot(data = transit.times[transit.times$o_microtype!=7,], aes(x = mode_chosen, y = mean.waittime)) +
#       geom_boxplot(aes(color = start_time_bin)) +
#       theme_bw() +
#       facet_wrap(.~o_microtype)
#   }
# }


#so for the time being, use wait and access time as a function of start.time and mode and origin microtypes
# if things are missing (actually only 2 in this case) use the mean value. 
transit.times = dat %>% 
  dplyr::group_by(start_time_bin,mode_chosen,o_microtype) %>%
  dplyr::summarise(mean.access.time = mean(access_time),
                   mean.waittime = mean(wait_time)) 

transit.times.mean = dat %>%
  dplyr::group_by(start_time_bin,mode_chosen) %>%
  dplyr::summarise(mean2.access.time = mean(access_time),
                   mean2.waittime = mean(wait_time)) 

#######################################################
### now generate data for non-chosen modes ###
dat2 = dat %>% 
  dplyr::select(trip_indx, sampno, perno, tdtrpnum, # tour_indx from tripchains dataset
                o_microtype, 
                d_microtype,
                o_geotype, 
                d_geotype,
                mode_chosen,         
                trip_purpose,
                access_time, wait_time, inv_time,
                start_time_bin,
                trpmiles) %>%
  unite('person_indx', sampno:perno, sep = '_',remove = F) 

# 896100 trips

# Added variables from the tripchain dataset
# tmp = dat2[,c('trip_indx','person_indx','trip_purpose','start_time_bin','o_microtype','d_microtype','o_geotype','d_geotype','trpmiles')] 
# altm.dat = NULL
# for(i in 1:length(unique(dat2$mode_chosen))){
#   altm.dat = rbind(altm.dat,tmp) 
# }
# altm.dat = altm.dat %>%
#   arrange(trip_indx) 
# altm.dat$mode_chosen = rep(c("walk", "hv",   "taxi", "bus" , "rail", "bike"),dim(dat2)[1]) 
# 
# altm.dat = altm.dat %>%
#   left_join(transit.times) %>% 
#   left_join(mode.speed) %>% 
#   left_join(mode.speed.mean)%>% 
#   left_join(dat2) %>% 
#   left_join(transit.times.mean) 
# 
# # Replacing NA values for the tripchain variables
# altm.dat = altm.dat %>%
#   group_by(trip_indx)%>%
#   fill(tour, tour_seg, tourtype, trpcnt, stops, tour_flg, trip_chain, .direction = "updown")%>%
#   dplyr::ungroup()

# code up the choices and generate travel times for alternative modes
# TravelMode = altm.dat %>%
#   dplyr::mutate(choice = case_when(
#     is.na(access_time) ~ FALSE,
#     TRUE ~ TRUE
#   )) %>%
#   dplyr::mutate(access_time = case_when(
#     is.na(access_time) ~ as.numeric(mean.access.time),
#     TRUE ~ as.numeric(access_time)
#   )) %>%
#   dplyr::mutate(access_time = case_when(
#     is.na(access_time) ~ as.numeric(mean2.access.time), # if by o_microtype not available, use average by mode and by starttime to impute
#     TRUE ~ as.numeric(access_time)
#   )) %>%
#   dplyr::mutate(wait_time = case_when(
#     is.na(wait_time) ~ as.numeric(mean.waittime),
#     TRUE ~ as.numeric(wait_time)
#   )) %>%
#   dplyr::mutate(wait_time = case_when(
#     is.na(wait_time) ~ as.numeric(mean2.waittime),
#     TRUE ~ as.numeric(wait_time)
#   )) %>%
#   dplyr::mutate(inv_time = case_when(
#     is.na(inv_time) & is.na(mean.speed) ~ as.numeric(trpmiles)/as.numeric(mean2.speed),
#     is.na(inv_time) & !(is.na(mean.speed)) ~ as.numeric(trpmiles)/as.numeric(mean.speed),
#     TRUE ~ as.numeric(inv_time)
#   )) %>%
#   dplyr::rename(mode = mode_chosen) %>%
#   dplyr::select(
#     -mean.access.time,
#     -mean.waittime,
#     -mean2.access.time,
#     -mean2.waittime,
#     -mean.speed,
#     -mean.speed.m_h,
#     - mean2.speed,
#     -mean2.speed.m_h,
#     -sampno,-perno, -tdtrpnum)

# merge to person weights
wts = wts %>%
  unite('person_indx', sampno:persno, sep = '_') 
TravelMode = dat2 %>%
  left_join(wts) # Joining, by = "person_indx"

# merge to dist_to_train 
# colnames(origin_train_dist) = tolower(colnames(origin_train_dist))
# #rename to CA convention
# origin_train_dist = origin_train_dist %>%
#   dplyr::rename(sampno = houseid,
#                 perno = personid) %>%
#   dplyr::select(-orig_country,-dest_country,-o_geotype,-o_microtype,-d_geotype,-d_microtype)%>% 
#   unite('trip_indx',sampno:tdtrpnum, sep="_") 
# 
# TravelMode = TravelMode %>%
#   left_join(origin_train_dist) # Joining, by = "trip_indx"

### prepare travel distance bins:
# within 1.3 walking distance
# 1.3 to 3 miles, biking distance
# 3 to 8 short driving
# 8 above long driving

# TravelMode = TravelMode %>%
#   dplyr::filter(!is.na(trpmiles))%>% 
#   dplyr::mutate(trip_dist_bin = case_when(
#     trpmiles <=1.3  ~ 's',
#     trpmiles >1.3 & trpmiles<=3 ~ 'm',
#     trpmiles>3 & trpmiles <=8 ~ 'l',
#     TRUE ~ 'll')) 


########################### merge and imput monetary cost ########################
tract <- read.csv(file =file.path(datadir, 'tripct.csv')) 
tract <- tract %>% 
  unite('trip_indx', HOUSEID:TDTRPNUM, sep = '_', remove = F)
tract$o_tract <- do.call(paste0, tract[c("ORIG_ST", "ORIG_CNTY", "ORIG_CT")]) 
tract$d_tract <- do.call(paste0, tract[c("DEST_ST", "DEST_CNTY", "DEST_CT")]) 
tract <- tract %>% 
  dplyr::select(trip_indx, o_tract, d_tract)

TravelMode <- TravelMode %>% 
  merge(tract, by="trip_indx") 
rm(tract)

dat <- TravelMode %>% 
  unite('o_microgeotype', c("o_geotype", "o_microtype"), sep = '_', remove = F) %>%
  dplyr::select(trip_indx, o_tract, o_microgeotype, o_microtype, o_geotype) %>% 
  distinct() 

####taxi fare
# taxi$tractcode <- as.character(taxi$tractcode)
dat$o_tract <- as.numeric(dat$o_tract)
taxi_dat <- taxi %>% merge(dat, by.x="tractcode", by.y="o_tract", all.y=TRUE)

taxi_dat_mean <- taxi_dat %>% 
  dplyr::group_by(o_microgeotype) %>% 
  dplyr::summarise(mean.minimumfare=mean(minimumfare, na.rm=TRUE),
                   mean.priceperminute=mean(priceperminute, na.rm=TRUE),
                   mean.pricepermile=mean(pricepermile, na.rm=TRUE))

taxi_dat_mean_m <- taxi_dat %>%
  dplyr::group_by(o_microtype) %>% 
  dplyr::summarise(mean.minimumfare.m=mean(minimumfare, na.rm=TRUE),
                   mean.priceperminute.m=mean(priceperminute, na.rm=TRUE),
                   mean.pricepermile.m=mean(pricepermile, na.rm=TRUE))

taxi_dat <- taxi_dat %>% 
  merge(taxi_dat_mean, by="o_microgeotype") %>% 
  merge(taxi_dat_mean_m, by="o_microtype") %>% 
  dplyr::mutate(minimumfare = case_when(is.na(minimumfare) ~ as.numeric(mean.minimumfare), TRUE ~ as.numeric(minimumfare)),
                priceperminute = case_when(is.na(priceperminute) ~ as.numeric(mean.priceperminute), TRUE ~ as.numeric(priceperminute)),
                pricepermile = case_when(is.na(pricepermile) ~ as.numeric(mean.pricepermile), TRUE ~ as.numeric(pricepermile))) %>%
  dplyr::mutate(minimumfare = case_when(is.na(minimumfare) ~ as.numeric(mean.minimumfare.m), TRUE ~ as.numeric(minimumfare)),
                priceperminute = case_when(is.na(priceperminute) ~ as.numeric(mean.priceperminute.m), TRUE ~ as.numeric(priceperminute)),
                pricepermile = case_when(is.na(pricepermile) ~ as.numeric(mean.pricepermile.m), TRUE ~ as.numeric(pricepermile))) %>%
  select(trip_indx, minimumfare, priceperminute, pricepermile) %>% mutate(mode="taxi")

rm(taxi, taxi_dat_mean, taxi_dat_mean_m)

# dat <- TravelMode %>% 
#   filter(choice==TRUE) %>% 
#   select(-mode)

taxi_dat <- taxi_dat %>% 
  merge(TravelMode, by="trip_indx") %>% 
  dplyr::mutate(taxi_cost=inv_time*priceperminute+trpmiles*pricepermile)

taxi_dat$taxi_cost[taxi_dat$taxi_cost<taxi_dat$minimumfare]<- taxi_dat$minimumfare[taxi_dat$taxi_cost<taxi_dat$minimumfare]
taxi_dat <- taxi_dat %>% 
  select(trip_indx, taxi_cost) %>% 
  dplyr::mutate(mode="taxi")

####public transit fare
dat <- TravelMode %>% 
  unite('o_microgeotype', c("o_geotype", "o_microtype"), sep = '_', remove = F) %>%
  select(trip_indx, mode_chosen, o_tract, o_microgeotype, o_microtype, o_geotype) %>% filter(mode_chosen=="rail" | mode_chosen=="bus")

# # CA only
# if(CA.TF){
# public <- public %>% 
#   dplyr::mutate(state=substr(tractcode,1,2)) %>% 
#   dplyr::filter(state=="06") %>% 
#   dplyr::select(-state)
# }

public$mode[public$mode=="rail_l"] <- "rail"
public$mode[public$mode=="rail_c"] <- "rail"
public <- public %>% 
  dplyr::group_by(mode, tractcode) %>% 
  dplyr::summarise(fare=mean(fare, na.rm=TRUE), mini=mean(mini, na.rm=TRUE), maxi=mean(maxi, na.rm=TRUE))

public_dat <- public %>% 
  merge(dat, by.x=c("tractcode", "mode"), by.y=c("o_tract", "mode_chosen"), all.y=TRUE)
public_dat_mean <- public_dat %>% 
  dplyr::group_by(mode, o_microgeotype) %>% 
  dplyr::summarise(mean.fare=mean(fare, na.rm=TRUE),
                   mean.mini=mean(mini, na.rm=TRUE),
                   mean.maxi=mean(maxi, na.rm=TRUE))
public_dat_mean_m <- public_dat %>%
  dplyr::group_by(mode, o_microtype) %>% 
  dplyr::summarise(mean.fare.m=mean(fare, na.rm=TRUE),
                   mean.mini.m=mean(mini, na.rm=TRUE),
                   mean.maxi.m=mean(maxi, na.rm=TRUE))
public_dat <- public_dat %>% 
  merge(public_dat_mean, by=c("o_microgeotype", "mode")) %>% 
  merge(public_dat_mean_m, by=c("o_microtype", "mode")) %>% 
  dplyr::mutate(fare = case_when(is.na(fare) ~ as.numeric(mean.fare), TRUE ~ as.numeric(fare)),
                mini = case_when(is.na(mini) ~ as.numeric(mean.mini), TRUE ~ as.numeric(mini)),
                maxi = case_when(is.na(maxi) ~ as.numeric(mean.maxi), TRUE ~ as.numeric(maxi))) %>%
  dplyr::mutate(fare = case_when(is.na(fare) ~ as.numeric(mean.fare.m), TRUE ~ as.numeric(fare)),
                mini = case_when(is.na(mini) ~ as.numeric(mean.mini.m), TRUE ~ as.numeric(mini)),
                maxi = case_when(is.na(maxi) ~ as.numeric(mean.maxi.m), TRUE ~ as.numeric(maxi))) %>%
  select(trip_indx, fare, mini, maxi, mode)

rm(public, public_dat_mean, public_dat_mean_m)

public_dat <- public_dat %>% 
  merge(TravelMode, by="trip_indx") %>%
  dplyr::mutate(public_cost = fare)

public_dat$public_cost[public_dat$mode=="rail" & is.na(public_dat$maxi)==FALSE & public_dat$trip_dist_bin=="ll"] <- public_dat$maxi[public_dat$mode=="rail" & is.na(public_dat$maxi)==FALSE & public_dat$trip_dist_bin=="ll"]
public_dat$public_cost[public_dat$mode=="rail" & is.na(public_dat$maxi)==FALSE & public_dat$trip_dist_bin=="l"] <- public_dat$maxi[public_dat$mode=="rail" & is.na(public_dat$maxi)==FALSE & public_dat$trip_dist_bin=="l"]

#adjust for senior fare 50% off
public_dat <- public_dat %>%
  select(trip_indx, public_cost, mode)

dat <- trips %>% 
  select(trip_indx, r_age_imp, trphhacc)

public_dat <- public_dat %>% 
  merge(dat, by="trip_indx")

public_dat$trphhacc[public_dat$trphhacc<0] <- 0
public_dat$trphhacc <- public_dat$trphhacc+1
public_dat$public_cost[public_dat$r_age_imp>=65] <- public_dat$public_cost[public_dat$r_age_imp>=65]*0.5+public_dat$public_cost[public_dat$r_age_imp>=65]*(public_dat$trphhacc[public_dat$r_age_imp>=65]-1)
public_dat <- public_dat %>% 
  select(trip_indx, public_cost, trphhacc, mode)

#####vehicle cost
vehicles <- vehicles %>% 
  unite('veh_indx', HOUSEID:VEHID, sep = '_', remove = F)

# prepare for fuel cost
trips <- trips %>% 
  unite('veh_indx', c("sampno", "vehid"), sep = '_', remove = F)

trip_fuel <- trips %>% 
  select(trip_indx, gasprice)

# process mpg data
trips_veh <- trips %>% 
  merge(vehicles, by="veh_indx", all.x=TRUE) %>% 
  dplyr::select(sampno, trip_indx, veh_indx, MAKE, MODEL, VEHTYPE)

veh_mpg <- read.csv(file=file.path(datadir, 'veh_mpg.csv'))

veh_mpg <- veh_mpg %>% 
  dplyr::select(makecode, modelcode, vehtype, mpg)
names(veh_mpg) <- c("MAKE", "MODEL", "VEHTYPE", "mpg")

trips_veh <- trips_veh %>% 
  left_join(veh_mpg, by=c("MAKE", "MODEL", "VEHTYPE"))

hhct <- read.csv(file=file.path(datadir, 'hhct.csv'))%>%
  mutate(HHCNTYFP = str_pad(HHCNTYFP, width =3, side = "left", pad = "0"),
         HHCT = str_pad(HHCT, width =6, side = "left", pad = "0"),
         GEOID = paste(HHSTFIPS, HHCNTYFP, HHCT, sep = ""),
         GEOID = str_pad(as.character(GEOID), width = 11, side = "left", pad = "0")) 

hhct <- hhct %>% 
  select(GEOID, HOUSEID)# %>%
# merge(pop , by = "GEOID") # merge population counts with household tracts
clust =  read.csv(file=file.path(datadir, "ccst_geoid_key_transp_geo_with_imputation.csv"))
clust <- clust %>%
  select(GEOID, microtype, geotype)
#merge trips by home tract
house_tract = hhct %>% 
  merge(clust, by = "GEOID") %>%
  mutate(h_geotype = geotype,
         h_microtype = microtype) %>% 
  unite('h_microgeotype', c("h_geotype", "h_microtype"), sep = '_', remove = F)

trips_veh <- trips_veh %>% 
  merge(house_tract, by.y="HOUSEID", by.x="sampno", all.x=T)

mpg_microgeo <- trips_veh %>% 
  dplyr::group_by(h_microgeotype) %>% 
  dplyr::summarise(mean.mpg=mean(mpg, na.rm=T))

trips_veh <- trips_veh %>% 
  merge(mpg_microgeo, by="h_microgeotype") %>%
  dplyr::select(trip_indx, mpg, mean.mpg) 


trips_veh$mpg[is.na(trips_veh$mpg)] <- trips_veh$mean.mpg[is.na(trips_veh$mpg)]

# parking cost
trips <- trips %>% 
  dplyr::mutate(parkingtime=case_when(strttime >=900 & strttime <1900 ~1, TRUE~0)) %>% 
  arrange(sampno, perno, tdaydate, strttime) %>%
  dplyr::group_by(sampno, perno, tdaydate) %>%
  dplyr::mutate(parking_dur = round((strttime - lag(strttime, default = strttime[1]))/100))


parking$state <- floor(parking$tractcode/(10^9))
parking$tractcode <- as.character(parking$tractcode)
parking$tractcode[parking$state<10] <- paste("0", parking$tractcode[parking$state<10], sep="")

dat <- TravelMode %>% 
  unite('d_microgeotype', c("d_geotype", "d_microtype"), sep = '_', remove = F) %>%
  select(trip_indx, d_tract, d_microgeotype, d_microtype, d_geotype) %>% distinct()

trips <- trips %>% 
  merge(dat, by.x="trip_indx")

# Note NP: would be better to rename some of these temp files so that the code isnt so senstive to the order
dat <- TravelMode %>% 
  select(trip_indx, trip_purpose) %>% 
  distinct()

trips <- trips %>% 
  merge(dat, by.x="trip_indx")

trips_parking <- trips %>% 
  dplyr::select(trip_indx, d_tract, parkingtime, parking_dur, trip_purpose, tdwknd) %>% 
  merge(parking, by.x=c("d_tract"), by.y=c("tractcode"), all.x=TRUE)

trips_parking$parking[is.na(trips_parking$parking)] <- 0

trips_parking <- trips_parking %>% 
  dplyr::select(trip_indx, parking, parkingtime, parking_dur, trip_purpose, tdwknd) %>% 
  dplyr::mutate(parking=parking*parking_dur) %>% 
  dplyr::mutate(parking=case_when(
    parking_dur==0 | tdwknd==1 | trip_purpose=="home" | trip_purpose=="work" ~ 0,
    TRUE~parking)) %>% 
  dplyr::select(trip_indx, parking)

# calculate the total vehicle costs
# dat <- TravelMode %>% 
#   filter(choice==TRUE)

hv_dat <- TravelMode %>% 
  merge(trips_parking, by="trip_indx") %>%
  merge(trips_veh, by="trip_indx") %>% 
  merge(trip_fuel, by="trip_indx") %>% 
  dplyr::mutate(veh_cost=gasprice*trpmiles/(100*mpg)+parking) %>% 
  dplyr::select(trip_indx, veh_cost) %>% 
  dplyr::mutate(mode="hv")

# # group by microtype
# hv_dat <- hv_dat %>% mutate(vehmile=gasprice/(100*mpg))
# hv_dat <- hv_dat %>% unite('o_microgeotype', c("o_geotype", "o_microtype"), sep = '_', remove = F)
# hv_dat <- hv_dat %>% unite('d_microgeotype', c("d_geotype", "d_microtype"), sep = '_', remove = F)
# hv_mile <- hv_dat %>% group_by(o_microgeotype) %>% summarise(hv_mile=mean(vehmile))
# hv_park <- hv_dat %>% group_by(d_microgeotype) %>% summarise(hv_park=mean(parking))


################merge all monetary costs
TravelMode <- TravelMode %>%
  dplyr::rename(mode = mode_chosen)
TravelMode <- TravelMode %>% 
  merge(hv_dat, by=c("trip_indx", "mode"), all.x=T) %>% 
  merge(public_dat, by=c("trip_indx", "mode"), all.x=T) %>% 
  merge(taxi_dat, by=c("trip_indx", "mode"), all.x=T)

TravelMode <- TravelMode %>% 
  dplyr::mutate(cost=case_when(mode=="hv"~veh_cost, mode=="bike"~0, mode=="walk"~0, mode=="bus"~public_cost,
                               mode=="rail"~public_cost, mode=="taxi"~taxi_cost))
names(TravelMode)
TravelMode <- TravelMode %>% 
  dplyr::select(-veh_cost, -public_cost, -trphhacc, -taxi_cost)

dat <- trips %>% 
  dplyr::select(trip_indx, trphhacc, hhvehcnt)

TravelMode <- TravelMode %>% 
  merge(dat, by="trip_indx")

TravelMode <- TravelMode %>% rename(HOUSEID = sampno)

location <- location %>% 
  dplyr::mutate(hometract = HHSTFIPS * 10^9 + HHCNTYFP*10^6 + HHCT) %>% 
  dplyr::select(HOUSEID, hometract)

# Update 10/25 adding bikeshare; 
TravelMode <- TravelMode %>% 
  merge(location, by="HOUSEID", all.x=T) %>% 
  merge(bikeshare, by.x="hometract", by.y="geoid", all.x=T)

# add home attributes
house_tract <- house_tract %>% 
  select(HOUSEID, h_geotype, h_microtype)

TravelMode <- TravelMode %>% 
  merge(house_tract, by="HOUSEID", all.x=T) 

TravelMode <- TravelMode %>%
  select(-hometract, -o_tract, -d_tract, -start_time_bin, -trphhacc, -density_land)

person_label <- read.csv(file=file.path(datadir, "nhts_user_classes_inc_veh_sr.csv"))
person_label <- person_label %>%
  select(person_indx, age, income, PopulationGroupID)

TravelMode <- TravelMode %>% 
  merge(person_label, by="person_indx", all.x=T) 

TravelMode <- TravelMode %>% rename(PERSONID = perno)

# save(TravelMode, file = file.path(datadir,'travelmode5_onerail_disability.RData'))

write.csv(TravelMode, file = file.path(datadir,'NHTS_data_with_time_cost.csv'))
# load(file = file.path(datadir,'travelmode5_onerail_disability_TpTr.RData'))

