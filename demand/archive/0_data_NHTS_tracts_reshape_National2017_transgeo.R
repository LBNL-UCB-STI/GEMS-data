## Rscript to process the trips data from National NHTS 2017
# Ling Jin 6/18/2020
# code up the mode, reshape the mode specific travel time variables.
# 6/24/2020 update: try estimate speed as a funciton of start.time, od-microtypes,and mode
#                   also try estimate access time and wait time as a function of start.time, origin-microtype, and mode
# 6/25/2020 update: split rail to 15-->commuter_rail; 16--> light_rail
# 6/29/2020 update: unobserved access and wait time imputed for rail modes for microtype5 using observed mean, 
#                   so that they are coded as a function of microtype, mode, and trip start time.
# 7/8/2020 update: use national data
# 7/15/2020 update: merge in person weights
# 7/20/2020 update: morning rush hour defined as starttime in [6am, 9am)
# 7/21/2020 update: prepare travel distance bins and interaction of constant with these bins
# 7/22/2020 update: only keep city bus trips for bus mode
# 3/15/2021 update: imputing based on transgeo results
# 5/20/2021 update: creating 24 hourly time bins instead of only 3
######################################################################
# mywd <- "/Users/lingjin/Dropbox/Research/FHWA_GeoType/Work/Rscripts"

mywd = "C:/FHWA/For FHWA folks/trip_generation_and_gems_inputs"

setwd(mywd)
source('Code/initialization.R')
source('Code/functions.R')
#figuredir <- "Figures"
outdir  <-  "CleanData"
#rdatadir <-  "../RData"
#datdir = '../Data/NHTS/nhts17_caltrans_tsdc_download/survey_data'
datdir = 'RawData'

library(stringr)
library(reshape2)
library(dplyr)
library(tidyr)
library(mlogit)
library(fastDummies)

#### local functions #########
unfold.list<- function(dat){
  dat$o_geotype = as.vector((do.call(rbind,lapply(dat$o_geotype,function(x)x[[1]]))))
  dat$d_geotype = as.vector((do.call(rbind,lapply(dat$d_geotype,function(x)x[[1]]))))
  dat$o_microtype = as.vector((do.call(rbind,lapply(dat$o_microtype,function(x)x[[1]]))))
  dat$d_microtype = as.vector((do.call(rbind,lapply(dat$d_microtype,function(x)x[[1]]))))
  dat = dat[complete.cases(dat),]
  return(dat)
}

nhts.lookup.var<- function(var, table = NULL, year = 2017){
  # function to look up values of a variable in nhts data
  # note that 2009 data has no excel codebook, so there is not searchable dictionary to use
  if(year == 2017){
    load('RawData/NHTS-Public/nhts.codebook.17.RData')
    lookuptable = lookuptable17
  }else{
    stop('no year found')
  }
  if(var %in% unique(lookuptable$name)){
    if(is.null(table)){
      return(lookuptable[lookuptable$name == var,])
    }else{
      if(table %in% unique(lookuptable$table)){
        return(lookuptable[lookuptable$name == var & lookuptable$table == table,])
      }else{
        stop('no table found')
      }
    }
  }else{
    stop('no variable found')
  }
}

map_labels <- function(dat,raw_colname = 'o_microtype',microtype = T){
  # function to map the raw cluster labels to relabeled names for micro and geotypes
  # input 
  #   dat is a dataframe that contains the raw cluster label in raw_colname
  #   if microtype = T, then take the raw label and map it to microtype labels, otherwise, map to geotype labels
  # return the updated label vector
  require(tidyr); require(dplyr)
  
  dat$inlabel = dat[,raw_colname]
  if(microtype){
    dat = dat %>%
      mutate(updatelabel = case_when(
        as.numeric(as.character(inlabel)) == 1 ~ as.numeric(4),
        as.numeric(as.character(inlabel)) == 2 ~ as.numeric(2),
        as.numeric(as.character(inlabel)) == 3 ~ as.numeric(6),
        as.numeric(as.character(inlabel)) == 4 ~ as.numeric(3),
        as.numeric(as.character(inlabel)) == 5 ~ as.numeric(5),
        as.numeric(as.character(inlabel)) == 6 ~ as.numeric(1),
        as.numeric(as.character(inlabel)) == 7 ~ as.numeric(0)
      ))
  }else{
    dat = dat %>%
      mutate(updatelabel = case_when(
        as.numeric(as.character(inlabel)) == 1 ~ as.character('G'),
        as.numeric(as.character(inlabel)) == 2 ~ as.character('E'),
        as.numeric(as.character(inlabel)) == 3 ~ as.character('H'),
        as.numeric(as.character(inlabel)) == 4 ~ as.character('D'),
        as.numeric(as.character(inlabel)) == 5 ~ as.character('B'),
        as.numeric(as.character(inlabel)) == 6 ~ as.character('F'),
        as.numeric(as.character(inlabel)) == 7 ~ as.character('C'),
        as.numeric(as.character(inlabel)) == 8 ~ as.character('I'),
        as.numeric(as.character(inlabel)) == 9 ~ as.character('A')
      ))
    
    
  }
  return(dat$updatelabel)
}


#####load the data and remove entries with NA ids ###############################

trips = read.csv(file = file.path(datdir,'NHTS-Public/trippub.csv'))
colnames(trips) = tolower(colnames(trips))
names(trips)
dim(trips) # 923572 trips

#rename to CA convention
trips = trips %>%
  rename(sampno = houseid,
         perno = personid)

trips = trips %>% 
  filter(!is.na(sampno) & !is.na(perno) & !is.na(tdtrpnum))
dim(trips) # same as before.

wts = read.csv(file = file.path(datdir,'NHTS-Public/perwgt.csv'))
# use WTPERFIN - person final weights
wts = wts %>%
  select(HOUSEID, PERSONID,  WTPERFIN) %>%
  rename(sampno = HOUSEID,
         persno = PERSONID,
         wtperfin = WTPERFIN)
##############################################
# 1. extract relevant variables
select.vars = 
  c("sampno", # household id
    "perno",  # person id within the household
    "tdtrpnum", # trip id
    "strttime", # trip start time
    'endtime', # trip end time
    "trvlcmin", # trip duration or travel time (calculated travel time)
    "trpmiles" , # total trip length in miles, or travel distance
    "trptrans", # derived trip mode, only single mode is selected.
    "trwaittm", # transit wait time in minutes, 0-60, set other values to NA
    "tracctm" , # trip time to transit station in minutes
    "drop_prk" , # park or dropped at transit station, 1-parked, 2-dropped, 3, bike/walked, other -NA
    "tregrtm" , # time used from transit to destination in minutes [0,180]
#    "tdwknd" , # is it a weekend trip. Note that generally there is only weekday data
#    "vmt_mile" ,# Trip distance in miles for personally driven vehicle trips, derived from route geometry returned by Google Maps API
    "whytrp1s" # trip purpose summary
) 

dat = trips[,select.vars]
  

#### take a look at wait_time by different bus modes ####
waittime = dat %>%
  filter(trwaittm>0)%>%
  group_by(trptrans) %>%
  summarise(mean_wait_time = mean(trwaittm))

######## only keep the city bus mode as bus mode #####
dat = dat %>%
  filter(!(trptrans %in% c(10,12,13,14))) %>%
  filter(!(trptrans == 11 & trwaittm >50))  # exclude bus trips with more than 20 min waiting time, 10% of the bus trips
  
######################################################
# 2. code up mode_chosen and upper_level_nest
dat = dat %>% 
  mutate(mode_chosen = case_when(
    trptrans %in% c(3,4,5,6,18,9) ~ 'hv', # human driving vehicle
    trptrans %in% c(2) ~ 'bike',
    trptrans %in% c(7,8) ~ 'scooter',
    trptrans %in% c(1) ~ 'walk',
    trptrans %in% 10:14 ~ 'bus', # now we only have 11 in the data
    trptrans %in% c(15) ~ 'rail_c',
    trptrans %in% c(16) ~ 'rail_l',
    trptrans %in% c(17) ~ 'taxi',
    TRUE ~ 'other'
  )) %>%
  mutate(upper_level_nest = case_when(
    mode_chosen %in% c('hv') ~ 'auto',
    mode_chosen %in% c('bike','scooter','walk') ~ 'micromobility',
    mode_chosen %in% c('bus','rail_l','rail_c','taxi') ~ 'transit',
    TRUE ~ 'other'
  ))


dat = as.data.frame(dat)

################################
# 3. code up trip purposes

for(var in c('whytrp1s')){
  tmp = nhts.lookup.var(var)
  if(class(dat[,var]) %in% c('integer','numeric')){
    dat[,var] = factor(dat[,var],levels = as.numeric(tmp$value), labels = tmp$label)
  }else{
    dat[,var] = factor(dat[,var],levels = tmp$value, labels = tmp$label)
  }
}

dat = dat %>%
  mutate(trip_purpose = case_when(
    whytrp1s %in% c('Social/Recreational') ~ 'social',
    whytrp1s %in% c('Home') ~ 'home',
    whytrp1s %in% c('Work') ~ 'work',
    whytrp1s %in% c('Shopping/Errands') ~ 'shopping',
    whytrp1s %in% c('Transport someone') ~ 'transp_someone',
    whytrp1s %in% c('Medical/Dental services') ~ 'medical',
    whytrp1s %in% c('Meals') ~ 'meals',
    whytrp1s %in% c('School/Daycare/Religious activity') ~ 'school',
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
  mutate(total_travel_time = as.numeric(trvlcmin)) %>% # total travel time in min
  mutate(access_time = case_when(  # travel time to transit station in minutes, if the mode is hv, bike/walk, this is 0
#    upper_level_nest %in% c('auto','micromobility') ~ 0,                
#    upper_level_nest == 'transit' &  tracctm == -1 ~ 0,
    upper_level_nest == 'transit' &  tracctm >0 ~ as.numeric(tracctm) ,
    TRUE ~ 0 
  )) %>%
  mutate(wait_time = case_when( ## wait time at transit stations in minutes, if the mode is hv, bike/walk, this is 0
    upper_level_nest == 'transit' &  trwaittm >0 ~ as.numeric(trwaittm) ,
    TRUE ~ 0
  )) %>%
  mutate(egress_time = case_when(
    upper_level_nest == 'transit' &  tregrtm >0 ~ as.numeric(tregrtm),
    TRUE ~ 0 
  ))%>%
  mutate(inv_time = total_travel_time - wait_time - access_time,
         inv_time2 = total_travel_time - wait_time - access_time - egress_time) %>% # in transit time if travel by transit or travel time if other modes
  filter(inv_time > 0) %>%  # remove transit trips with derived in transit time less than 0     
  filter(trpmiles >0) %>% # excluding 91 trips that had NA in travel distance
  mutate(inv_speed = as.numeric(trpmiles)/inv_time) %>%
  filter(inv_speed < 100/60) # removed 257 trips with derived speed exceedingly large (i.e. >100 mile/h)

###################################################
# 5. Generate travel time bins (now not decided): morning rush hour [6am, 9am), evening rush hour [4pm, 7pm], and other hours
# NOTE NDP 5/20/21: We need to preserve the variable "strttime" so we can create time bins on the fly
dat = dat %>%
  mutate(start_time_bin = case_when(
    strttime >=600 & strttime <900 ~ 'morning_rush',
    strttime >=1600 & strttime <=1900 ~ 'evening_rush',
    TRUE ~ 'other_time'
  ))




#######################################################################
# 6. to generate variables for alternative modes, need derive travel speed by each mode, then summarize by micro_geotype o-d pairs and by travel start-time bin
#    then generate travel time or inv_time for each mode
master.dat = dat

load(file = file.path(outdir,'NHTS/NHTS_tract_od.tr.purpose2.National_transgeo.RData')) # load the trip purpose data

#od.tr.purpose2 = unfold.list(od.tr.purpose2)
# map to our final cluster labels (transgeo version already mapped to final labels)
#od.tr.purpose2$o_microtype = map_labels(od.tr.purpose2,raw_colname = 'o_microtype',microtype = T)
#od.tr.purpose2$d_microtype = map_labels(od.tr.purpose2,raw_colname = 'd_microtype',microtype = T)
#od.tr.purpose2$o_geotype = map_labels(od.tr.purpose2,raw_colname = 'o_geotype',microtype = F)
#od.tr.purpose2$d_geotype = map_labels(od.tr.purpose2,raw_colname = 'd_geotype',microtype = F)

tmp = od.tr.purpose2 %>%
  unite('o_geo_mic', o_geotype, o_microtype, remove = F) %>%
  unite('d_geo_mic', d_geotype, d_microtype, remove =F) %>%
  unite('od_geo_mic', o_geo_mic, d_geo_mic,remove = F) %>%
  rename(sampno = HOUSEID,
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
testing = F
if(testing){
  mode.speed = dat %>%
    group_by(start_time_bin,od_geo_mic, mode_chosen,o_microtype,d_microtype) %>%
    summarise(mean.speed = mean(inv_speed),
              n.sample = n())
  
  mode.speed %>%
    group_by(start_time_bin, mode_chosen) %>%
    summarise(N_unique_od_pairs = n_distinct(od_geo_mic))
  
  
  mode.speed %>%
    group_by( mode_chosen) %>%
    summarise(N_unique_od_pairs = n_distinct(od_geo_mic))
  
  # how about group by od-micriotype pairs --> We can see there are a number of locations with no taxi and rail speed observed
  groupby.od.micro = F
  if(groupby.od.micro ){
    mode.speed = dat %>%
    group_by(start_time_bin,o_microtype,d_microtype, mode_chosen) %>%
    summarise(mean.speed = mean(inv_speed),
              n.sample = n())
  
  ggplot(data = mode.speed, aes(x = mode_chosen, y = mean.speed)) +
    geom_boxplot(aes(color = start_time_bin)) +
    theme_bw() +
    facet_grid(o_microtype ~ d_microtype)
  }
  
  
  # look at the mode speed (first averaged over micro-geotype pairs)
  ggplot(data = mode.speed, aes(x = mode_chosen, y = mean.speed*60)) +
    geom_boxplot(aes(color = start_time_bin))
  
  # look at the average raw speed 
  ggplot(data = dat, aes(x = mode_chosen, y = inv_speed*60)) +
    geom_boxplot(aes(color = start_time_bin))+
    theme_bw() +
    facet_grid(o_microtype ~ d_microtype)
}


# So for the time being, generate average speed by mode 
# as a function of start.time and mode, and od microtype pairs
# if no observation at od microptype pairs, use the mean value.
# LJ 8/12/2020 change to harmonic means harmonic.mean(x)

mode.speed = dat %>%
  group_by(start_time_bin, mode_chosen,o_microtype,d_microtype) %>%
  summarise(mean.speed = harmonic.mean(inv_speed),
            mean.speed.m_h = harmonic.mean(inv_speed) * 60)

mode.speed.mean = dat %>%
  group_by(start_time_bin, mode_chosen) %>%
  summarise(mean2.speed = harmonic.mean(inv_speed),
            mean2.speed.m_h = harmonic.mean(inv_speed) * 60 )

# observed access.time and wait.time by origin, start.time. and mode chosen
# full set should include 7*9 (origin micro-geo types) * 3 (time bin) * 6 (modes) = 1134
# instead we observe 43 origin micro-geo types
# there are more origin types in hv than other modes, so observed wait time and access time will not cover all the origins in the auto (hv) trips
#so for the time being, use wait and access time as a function of start.time and mode. 
testing = F
if(testing){
  transit.times = dat %>% 
    group_by(start_time_bin,o_geo_mic, mode_chosen) %>%
    summarise(mean.access.time = mean(access_time),
              mean.waittime = mean(wait_time))
  
  transit.times %>%
    group_by(start_time_bin,mode_chosen) %>%
    summarise(N_origin = n_distinct(o_geo_mic)) 
  
  transit.times %>%
    group_by(mode_chosen) %>%
    summarise(N_origin = n_distinct(o_geo_mic))
  
  ggplot(data = transit.times, aes(x = mode_chosen, y = mean.access.time)) +
    geom_boxplot(aes(color = start_time_bin))
  
  ggplot(data = transit.times, aes(x = mode_chosen, y = mean.waittime)) +
    geom_boxplot(aes(color = start_time_bin))
  
  # how about only as a function of origin microtype rather than micro-geotype combinations?
  # --> we now have everything observed. 
  
  groupby.o.micro = F
  if(groupby.o.micro){
    transit.times = dat %>% 
      group_by(start_time_bin,o_microtype, mode_chosen) %>%
      summarise(mean.access.time = mean(access_time),
                mean.waittime = mean(wait_time))
    
    ggplot(data = transit.times[transit.times$o_microtype!=7,], aes(x = mode_chosen, y = mean.access.time)) +
      geom_boxplot(aes(color = start_time_bin)) +
      theme_bw() +
      facet_wrap(.~o_microtype)
    
    ggplot(data = transit.times[transit.times$o_microtype!=7,], aes(x = mode_chosen, y = mean.waittime)) +
      geom_boxplot(aes(color = start_time_bin)) +
      theme_bw() +
      facet_wrap(.~o_microtype)
  }
}


#so for the time being, use wait and access time as a function of start.time and mode and origin microtypes
# if things are missing (actually only 2 in this case) use the mean value. 
transit.times = dat %>% 
  group_by(start_time_bin,mode_chosen,o_microtype) %>%
  summarise(mean.access.time = mean(access_time),
            mean.waittime = mean(wait_time))

transit.times.mean = dat %>% 
  group_by(start_time_bin,mode_chosen) %>%
  summarise(mean2.access.time = mean(access_time),
            mean2.waittime = mean(wait_time))

#######################################################
### now generate data for non-chosen modes ###
dat2 = dat %>% 
  select(sampno, perno, tdtrpnum,
         o_microtype, 
         d_microtype,
         o_geotype, 
         d_geotype,
         mode_chosen,         
         trip_purpose,
         access_time, wait_time, inv_time, 
         start_time_bin,
         strttime,
         trpmiles) %>%
  unite('trip_indx', sampno:tdtrpnum, sep = '_', remove = F) %>%
  unite('person_indx', sampno:perno, sep = '_',remove = F)

tmp = dat2[,c('trip_indx','person_indx','trip_purpose','start_time_bin','o_microtype','d_microtype','o_geotype','d_geotype', 'strttime', 'trpmiles')] 
altm.dat = NULL
for(i in 1:length(unique(dat2$mode_chosen))){
  altm.dat = rbind(altm.dat,tmp)
}
altm.dat = altm.dat %>%
  arrange(trip_indx)
altm.dat$mode_chosen = rep(c("walk", "hv", "taxi", "bus" , "rail_l","rail_c", "bike"),dim(dat2)[1])

altm.dat = altm.dat %>%
  left_join(transit.times) %>%
  left_join(mode.speed) %>%
  left_join(mode.speed.mean)%>%
  left_join(dat2) %>%
  left_join(transit.times.mean)

# code up the choices and generate travel times for alternative modes
TravelMode = altm.dat %>%
  mutate(choice = case_when(
    is.na(access_time) ~ FALSE,
    TRUE ~ TRUE
  )) %>%
  mutate(access_time = case_when(
    is.na(access_time) ~ as.numeric(mean.access.time),
    TRUE ~ as.numeric(access_time)
  )) %>%
  mutate(access_time = case_when(
    is.na(access_time) ~ as.numeric(mean2.access.time), # if by o_microtype not available, use average by mode and by starttime to impute
    TRUE ~ as.numeric(access_time)
  )) %>%
  mutate(wait_time = case_when(
    is.na(wait_time) ~ as.numeric(mean.waittime),
    TRUE ~ as.numeric(wait_time)
  )) %>%
  mutate(wait_time = case_when(
    is.na(wait_time) ~ as.numeric(mean2.waittime),
    TRUE ~ as.numeric(wait_time)
  )) %>%
  mutate(inv_time = case_when(
    is.na(inv_time) & is.na(mean.speed) ~ as.numeric(trpmiles)/as.numeric(mean2.speed),
    is.na(inv_time) & !(is.na(mean.speed)) ~ as.numeric(trpmiles)/as.numeric(mean.speed),
    TRUE ~ as.numeric(inv_time)
  )) %>%
  rename(mode = mode_chosen) %>%
  select(
         -mean.access.time,
         -mean.waittime,
         -mean2.access.time,
         -mean2.waittime,
         -mean.speed,
         -mean.speed.m_h,
         - mean2.speed,
         -mean2.speed.m_h,
         -sampno,-perno, -tdtrpnum)


# merge to person weights
wts = wts %>%
  unite('person_indx', sampno:persno, sep = '_') 
TravelMode = TravelMode %>%
  left_join(wts)

save(TravelMode, file = file.path(outdir,'NHTS/TravelMode.tmp.RData'))
load(file = file.path(outdir,'NHTS/TravelMode.tmp.RData'))
# merge to dist_to_train 
# load the train distance data merged to trip locations.
load(file= file.path(outdir,'NHTS/NHTS_tract_origin_train_dist_transgeo.RData'))
colnames(origin_train_dist) = tolower(colnames(origin_train_dist))
#rename to CA convention
origin_train_dist = origin_train_dist %>%
  rename(sampno = houseid,
         perno = personid) %>%
  select(-orig_country,-dest_country,-o_geotype,-o_microtype,-d_geotype,-d_microtype)%>%
  unite('trip_indx',sampno:tdtrpnum, sep="_")

TravelMode = TravelMode %>%
  left_join(origin_train_dist)

### prepare travel distance bins:
# within 1.3 walking distance
# 1.3 to 3 miles, biking distance
# 3 to 8 short driving
# 8 above long driving

TravelMode = TravelMode %>%
  filter(!is.na(trpmiles))%>%
  mutate(trip_dist_bin = case_when(
    trpmiles <=1.3  ~ 's',
    trpmiles >1.3 & trpmiles<=3 ~ 'm',
    trpmiles>3 & trpmiles <=8 ~ 'l',
    TRUE ~ 'll'))

save(TravelMode, file = file.path(outdir,'NHTS/travelmode.test.data.national_transgeo.RData'))

load(file = file.path(outdir,'NHTS/travelmode.test.data.national_transgeo.RData'))

# if on databucket, prepare a version of data with trip start tract id. 
merge.tract = F
if(merge.tract){
  # nhts.dir <- 'C:/NHTS/2017_tripct/'
  
  triploc = read.csv(file = file.path(datdir, 'NHTS-Tract/tripct.csv') ) %>%
    unite('o_geoid',ORIG_ST : ORIG_CT, sep='') %>%
    unite('d_geoid',DEST_ST : DEST_CT, sep = '') %>%
    mutate(o_geoid = as.numeric(o_geoid),
           d_geoid = as.numeric(d_geoid))%>%
    unite(trip_indx, HOUSEID:TDTRPNUM, sep = '_',remove = F)%>%
    select(-ORIG_COUNTRY, - DEST_COUNTRY) %>%
    right_join(TravelMode)
  
TravelMode =  triploc
save(TravelMode, file = file.path(outdir,'NHTS/travelmode.data.national.withgeoid_transgeo.RData'))
#load(file = file.path('C:/Data/NHTS','travelmode.data.national.withgeoid.RData'))

}

###################################################################
# 7. As we have different specificatoin for each nest (transit is different from others), 
# we have to manually generate some variable and specify the model formula
###### We need to generate indicator variables for mode and trip purposes (later we will need pop classes) 
tmp = TravelMode %>%
  rename(trip_purp = trip_purpose)%>%
  fastDummies::dummy_cols(select_columns = c('trip_purp','mode','trip_dist_bin')) %>%
  mutate(# generate covariates for auto nest modes (current only one)
    
    # mode specific constant
    mode_hv = mode_hv,
    
    # mode specific constant interacting with distance bins
    mode_hvXtrip_s = mode_hv*trip_dist_bin_s,
    mode_hvXtrip_m = mode_hv*trip_dist_bin_m,
    mode_hvXtrip_l = mode_hv*trip_dist_bin_l,
    mode_hvXtrip_ll = mode_hv*trip_dist_bin_ll,
    
    # mode specific travel time
    mode_hvXinv_time = mode_hv*inv_time,
    
    #trip purpose and mode specific constant
    mode_hvXtrip_purp_home = mode_hv*trip_purp_home,
    mode_hvXtrip_purp_meals = mode_hv*trip_purp_meals,
    mode_hvXtrip_purp_medical = mode_hv*trip_purp_medical,
    mode_hvXtrip_purp_other = mode_hv*trip_purp_other,
    mode_hvXtrip_purp_school = mode_hv*trip_purp_school,
    mode_hvXtrip_purp_shopping = mode_hv*trip_purp_shopping,
    mode_hvXtrip_purp_social = mode_hv*trip_purp_social,
    mode_hvXtrip_purp_transp_someone = mode_hv*trip_purp_transp_someone,
    mode_hvXtrip_purp_work = mode_hv*trip_purp_work,
    
    # trip purpose and mode specific travel time
    mode_hvXtrip_purp_homeXinv_time = mode_hv*trip_purp_home*inv_time,
    mode_hvXtrip_purp_mealsXinv_time = mode_hv*trip_purp_meals*inv_time,
    mode_hvXtrip_purp_medicalXinv_time = mode_hv*trip_purp_medical*inv_time,
    mode_hvXtrip_purp_otherXinv_time = mode_hv*trip_purp_other*inv_time,
    mode_hvXtrip_purp_schoolXinv_time = mode_hv*trip_purp_school*inv_time,
    mode_hvXtrip_purp_shoppingXinv_time = mode_hv*trip_purp_shopping*inv_time,
    mode_hvXtrip_purp_socialXinv_time = mode_hv*trip_purp_social*inv_time,
    mode_hvXtrip_purp_transp_someoneXinv_time = mode_hv*trip_purp_transp_someone*inv_time,
    mode_hvXtrip_purp_workXinv_time = mode_hv*trip_purp_work*inv_time  ) %>%
  mutate( # generate walk mode variables
    # mode specific constant
    mode_walk = mode_walk,
    
    # mode specific constant interacting with distance bins
    mode_walkXtrip_s = mode_walk*trip_dist_bin_s,
    mode_walkXtrip_m = mode_walk*trip_dist_bin_m,
    mode_walkXtrip_l = mode_walk*trip_dist_bin_l,
    mode_walkXtrip_ll = mode_walk*trip_dist_bin_ll,
    
    # mode specific travel time
    mode_walkXinv_time = mode_walk*inv_time,
    
    #trip purpose and mode specific constant
    mode_walkXtrip_purp_home = mode_walk*trip_purp_home,
    mode_walkXtrip_purp_meals = mode_walk*trip_purp_meals,
    mode_walkXtrip_purp_medical = mode_walk*trip_purp_medical,
    mode_walkXtrip_purp_other = mode_walk*trip_purp_other,
    mode_walkXtrip_purp_school = mode_walk*trip_purp_school,
    mode_walkXtrip_purp_shopping = mode_walk*trip_purp_shopping,
    mode_walkXtrip_purp_social = mode_walk*trip_purp_social,
    mode_walkXtrip_purp_transp_someone = mode_walk*trip_purp_transp_someone,
    mode_walkXtrip_purp_work = mode_walk*trip_purp_work,
    
    # trip purpose and mode specific travel time
    mode_walkXtrip_purp_homeXinv_time = mode_walk*trip_purp_home*inv_time,
    mode_walkXtrip_purp_mealsXinv_time = mode_walk*trip_purp_meals*inv_time,
    mode_walkXtrip_purp_medicalXinv_time = mode_walk*trip_purp_medical*inv_time,
    mode_walkXtrip_purp_otherXinv_time = mode_walk*trip_purp_other*inv_time,
    mode_walkXtrip_purp_schoolXinv_time = mode_walk*trip_purp_school*inv_time,
    mode_walkXtrip_purp_shoppingXinv_time = mode_walk*trip_purp_shopping*inv_time,
    mode_walkXtrip_purp_socialXinv_time = mode_walk*trip_purp_social*inv_time,
    mode_walkXtrip_purp_transp_someoneXinv_time = mode_walk*trip_purp_transp_someone*inv_time,
    mode_walkXtrip_purp_workXinv_time = mode_walk*trip_purp_work*inv_time  ) %>%
  mutate( # generate bike related variables
    # mode specific constant
    mode_bike = mode_bike,

    # mode specific constant interacting with distance bins
    mode_bikeXtrip_s = mode_bike*trip_dist_bin_s,
    mode_bikeXtrip_m = mode_bike*trip_dist_bin_m,
    mode_bikeXtrip_l = mode_bike*trip_dist_bin_l,
    mode_bikeXtrip_ll = mode_bike*trip_dist_bin_ll,
    
    # mode specific travel time
    mode_bikeXinv_time = mode_bike*inv_time,
    
    #trip purpose and mode specific constant
    mode_bikeXtrip_purp_home = mode_bike*trip_purp_home,
    mode_bikeXtrip_purp_meals = mode_bike*trip_purp_meals,
    mode_bikeXtrip_purp_medical = mode_bike*trip_purp_medical,
    mode_bikeXtrip_purp_other = mode_bike*trip_purp_other,
    mode_bikeXtrip_purp_school = mode_bike*trip_purp_school,
    mode_bikeXtrip_purp_shopping = mode_bike*trip_purp_shopping,
    mode_bikeXtrip_purp_social = mode_bike*trip_purp_social,
    mode_bikeXtrip_purp_transp_someone = mode_bike*trip_purp_transp_someone,
    mode_bikeXtrip_purp_work = mode_bike*trip_purp_work,
    
    # trip purpose and mode specific travel time
    mode_bikeXtrip_purp_homeXinv_time = mode_bike*trip_purp_home*inv_time,
    mode_bikeXtrip_purp_mealsXinv_time = mode_bike*trip_purp_meals*inv_time,
    mode_bikeXtrip_purp_medicalXinv_time = mode_bike*trip_purp_medical*inv_time,
    mode_bikeXtrip_purp_otherXinv_time = mode_bike*trip_purp_other*inv_time,
    mode_bikeXtrip_purp_schoolXinv_time = mode_bike*trip_purp_school*inv_time,
    mode_bikeXtrip_purp_shoppingXinv_time = mode_bike*trip_purp_shopping*inv_time,
    mode_bikeXtrip_purp_socialXinv_time = mode_bike*trip_purp_social*inv_time,
    mode_bikeXtrip_purp_transp_someoneXinv_time = mode_bike*trip_purp_transp_someone*inv_time,
    mode_bikeXtrip_purp_workXinv_time = mode_bike*trip_purp_work*inv_time  ) %>%
  mutate( # generate transit-bus related variables
    # mode specific constant
    mode_bus = mode_bus,
    
    mode_busXtrip_s = mode_bus*trip_dist_bin_s,
    mode_busXtrip_m = mode_bus*trip_dist_bin_m,
    mode_busXtrip_l = mode_bus*trip_dist_bin_l,
    mode_busXtrip_ll = mode_bus*trip_dist_bin_ll,
    
    
    # mode specific travel time
    mode_busXinv_time = mode_bus*inv_time,
    mode_busXaccess_time = mode_bus*access_time,
    mode_busXwait_time = mode_bus*wait_time,
    
    
    #trip purpose and mode specific constant
    mode_busXtrip_purp_home = mode_bus*trip_purp_home,
    mode_busXtrip_purp_meals = mode_bus*trip_purp_meals,
    mode_busXtrip_purp_medical = mode_bus*trip_purp_medical,
    mode_busXtrip_purp_other = mode_bus*trip_purp_other,
    mode_busXtrip_purp_school = mode_bus*trip_purp_school,
    mode_busXtrip_purp_shopping = mode_bus*trip_purp_shopping,
    mode_busXtrip_purp_social = mode_bus*trip_purp_social,
    mode_busXtrip_purp_transp_someone = mode_bus*trip_purp_transp_someone,
    mode_busXtrip_purp_work = mode_bus*trip_purp_work,
    
    # trip purpose and mode specific inv travel time
    mode_busXtrip_purp_homeXinv_time = mode_bus*trip_purp_home*inv_time,
    mode_busXtrip_purp_mealsXinv_time = mode_bus*trip_purp_meals*inv_time,
    mode_busXtrip_purp_medicalXinv_time = mode_bus*trip_purp_medical*inv_time,
    mode_busXtrip_purp_otherXinv_time = mode_bus*trip_purp_other*inv_time,
    mode_busXtrip_purp_schoolXinv_time = mode_bus*trip_purp_school*inv_time,
    mode_busXtrip_purp_shoppingXinv_time = mode_bus*trip_purp_shopping*inv_time,
    mode_busXtrip_purp_socialXinv_time = mode_bus*trip_purp_social*inv_time,
    mode_busXtrip_purp_transp_someoneXinv_time = mode_bus*trip_purp_transp_someone*inv_time,
    mode_busXtrip_purp_workXinv_time = mode_bus*trip_purp_work*inv_time,
    
    # trip purpose and mode specific access  time
    mode_busXtrip_purp_homeXaccess_time = mode_bus*trip_purp_home*access_time,
    mode_busXtrip_purp_mealsXaccess_time = mode_bus*trip_purp_meals*access_time,
    mode_busXtrip_purp_medicalXaccess_time = mode_bus*trip_purp_medical*access_time,
    mode_busXtrip_purp_otherXaccess_time = mode_bus*trip_purp_other*access_time,
    mode_busXtrip_purp_schoolXaccess_time = mode_bus*trip_purp_school*access_time,
    mode_busXtrip_purp_shoppingXaccess_time = mode_bus*trip_purp_shopping*access_time,
    mode_busXtrip_purp_socialXaccess_time = mode_bus*trip_purp_social*access_time,
    mode_busXtrip_purp_transp_someoneXaccess_time = mode_bus*trip_purp_transp_someone*access_time,
    mode_busXtrip_purp_workXaccess_time = mode_bus*trip_purp_work*access_time,
    
    # trip purpose and mode specific wait  time
    mode_busXtrip_purp_homeXwait_time = mode_bus*trip_purp_home*wait_time,
    mode_busXtrip_purp_mealsXwait_time = mode_bus*trip_purp_meals*wait_time,
    mode_busXtrip_purp_medicalXwait_time = mode_bus*trip_purp_medical*wait_time,
    mode_busXtrip_purp_otherXwait_time = mode_bus*trip_purp_other*wait_time,
    mode_busXtrip_purp_schoolXwait_time = mode_bus*trip_purp_school*wait_time,
    mode_busXtrip_purp_shoppingXwait_time = mode_bus*trip_purp_shopping*wait_time,
    mode_busXtrip_purp_socialXwait_time = mode_bus*trip_purp_social*wait_time,
    mode_busXtrip_purp_transp_someoneXwait_time = mode_bus*trip_purp_transp_someone*wait_time,
    mode_busXtrip_purp_workXwait_time = mode_bus*trip_purp_work*wait_time) %>%
  mutate( # generate transit-taxi related variables
    # mode specific constant
    mode_taxi = mode_taxi,
    
    mode_taxiXtrip_s = mode_taxi*trip_dist_bin_s,
    mode_taxiXtrip_m = mode_taxi*trip_dist_bin_m,
    mode_taxiXtrip_l = mode_taxi*trip_dist_bin_l,
    mode_taxiXtrip_ll = mode_taxi*trip_dist_bin_ll,
    
    
    # mode specific travel time
    mode_taxiXinv_time = mode_taxi*inv_time,
    mode_taxiXaccess_time = mode_taxi*access_time,
    mode_taxiXwait_time = mode_taxi*wait_time,
    
    #trip purpose and mode specific constant
    mode_taxiXtrip_purp_home = mode_taxi*trip_purp_home,
    mode_taxiXtrip_purp_meals = mode_taxi*trip_purp_meals,
    mode_taxiXtrip_purp_medical = mode_taxi*trip_purp_medical,
    mode_taxiXtrip_purp_other = mode_taxi*trip_purp_other,
    mode_taxiXtrip_purp_school = mode_taxi*trip_purp_school,
    mode_taxiXtrip_purp_shopping = mode_taxi*trip_purp_shopping,
    mode_taxiXtrip_purp_social = mode_taxi*trip_purp_social,
    mode_taxiXtrip_purp_transp_someone = mode_taxi*trip_purp_transp_someone,
    mode_taxiXtrip_purp_work = mode_taxi*trip_purp_work,
    
    # trip purpose and mode specific inv travel time
    mode_taxiXtrip_purp_homeXinv_time = mode_taxi*trip_purp_home*inv_time,
    mode_taxiXtrip_purp_mealsXinv_time = mode_taxi*trip_purp_meals*inv_time,
    mode_taxiXtrip_purp_medicalXinv_time = mode_taxi*trip_purp_medical*inv_time,
    mode_taxiXtrip_purp_otherXinv_time = mode_taxi*trip_purp_other*inv_time,
    mode_taxiXtrip_purp_schoolXinv_time = mode_taxi*trip_purp_school*inv_time,
    mode_taxiXtrip_purp_shoppingXinv_time = mode_taxi*trip_purp_shopping*inv_time,
    mode_taxiXtrip_purp_socialXinv_time = mode_taxi*trip_purp_social*inv_time,
    mode_taxiXtrip_purp_transp_someoneXinv_time = mode_taxi*trip_purp_transp_someone*inv_time,
    mode_taxiXtrip_purp_workXinv_time = mode_taxi*trip_purp_work*inv_time,
    
    # trip purpose and mode specific access  time
    mode_taxiXtrip_purp_homeXaccess_time = mode_taxi*trip_purp_home*access_time,
    mode_taxiXtrip_purp_mealsXaccess_time = mode_taxi*trip_purp_meals*access_time,
    mode_taxiXtrip_purp_medicalXaccess_time = mode_taxi*trip_purp_medical*access_time,
    mode_taxiXtrip_purp_otherXaccess_time = mode_taxi*trip_purp_other*access_time,
    mode_taxiXtrip_purp_schoolXaccess_time = mode_taxi*trip_purp_school*access_time,
    mode_taxiXtrip_purp_shoppingXaccess_time = mode_taxi*trip_purp_shopping*access_time,
    mode_taxiXtrip_purp_socialXaccess_time = mode_taxi*trip_purp_social*access_time,
    mode_taxiXtrip_purp_transp_someoneXaccess_time = mode_taxi*trip_purp_transp_someone*access_time,
    mode_taxiXtrip_purp_workXaccess_time = mode_taxi*trip_purp_work*access_time,
    
    # trip purpose and mode specific wait  time
    mode_taxiXtrip_purp_homeXwait_time = mode_taxi*trip_purp_home*wait_time,
    mode_taxiXtrip_purp_mealsXwait_time = mode_taxi*trip_purp_meals*wait_time,
    mode_taxiXtrip_purp_medicalXwait_time = mode_taxi*trip_purp_medical*wait_time,
    mode_taxiXtrip_purp_otherXwait_time = mode_taxi*trip_purp_other*wait_time,
    mode_taxiXtrip_purp_schoolXwait_time = mode_taxi*trip_purp_school*wait_time,
    mode_taxiXtrip_purp_shoppingXwait_time = mode_taxi*trip_purp_shopping*wait_time,
    mode_taxiXtrip_purp_socialXwait_time = mode_taxi*trip_purp_social*wait_time,
    mode_taxiXtrip_purp_transp_someoneXwait_time = mode_taxi*trip_purp_transp_someone*wait_time,
    mode_taxiXtrip_purp_workXwait_time = mode_taxi*trip_purp_work*wait_time) %>%
  mutate( # generate transit-rail_l related variables
    # mode specific constant
    mode_rail_l = mode_rail_l,
    
    mode_rail_lXtrip_s = mode_rail_l*trip_dist_bin_s,
    mode_rail_lXtrip_m = mode_rail_l*trip_dist_bin_m,
    mode_rail_lXtrip_l = mode_rail_l*trip_dist_bin_l,
    mode_rail_lXtrip_ll = mode_rail_l*trip_dist_bin_ll,
    
    
    
    # mode specific travel time
    mode_rail_lXinv_time = mode_rail_l*inv_time,
    mode_rail_lXaccess_time = mode_rail_l*access_time,
    mode_rail_lXwait_time = mode_rail_l*wait_time,
    
    #trip purpose and mode specific constant
    mode_rail_lXtrip_purp_home = mode_rail_l*trip_purp_home,
    mode_rail_lXtrip_purp_meals = mode_rail_l*trip_purp_meals,
    mode_rail_lXtrip_purp_medical = mode_rail_l*trip_purp_medical,
    mode_rail_lXtrip_purp_other = mode_rail_l*trip_purp_other,
    mode_rail_lXtrip_purp_school = mode_rail_l*trip_purp_school,
    mode_rail_lXtrip_purp_shopping = mode_rail_l*trip_purp_shopping,
    mode_rail_lXtrip_purp_social = mode_rail_l*trip_purp_social,
    mode_rail_lXtrip_purp_transp_someone = mode_rail_l*trip_purp_transp_someone,
    mode_rail_lXtrip_purp_work = mode_rail_l*trip_purp_work,
    
    # trip purpose and mode specific inv travel time
    mode_rail_lXtrip_purp_homeXinv_time = mode_rail_l*trip_purp_home*inv_time,
    mode_rail_lXtrip_purp_mealsXinv_time = mode_rail_l*trip_purp_meals*inv_time,
    mode_rail_lXtrip_purp_medicalXinv_time = mode_rail_l*trip_purp_medical*inv_time,
    mode_rail_lXtrip_purp_otherXinv_time = mode_rail_l*trip_purp_other*inv_time,
    mode_rail_lXtrip_purp_schoolXinv_time = mode_rail_l*trip_purp_school*inv_time,
    mode_rail_lXtrip_purp_shoppingXinv_time = mode_rail_l*trip_purp_shopping*inv_time,
    mode_rail_lXtrip_purp_socialXinv_time = mode_rail_l*trip_purp_social*inv_time,
    mode_rail_lXtrip_purp_transp_someoneXinv_time = mode_rail_l*trip_purp_transp_someone*inv_time,
    mode_rail_lXtrip_purp_workXinv_time = mode_rail_l*trip_purp_work*inv_time,
    
    # trip purpose and mode specific access  time
    mode_rail_lXtrip_purp_homeXaccess_time = mode_rail_l*trip_purp_home*access_time,
    mode_rail_lXtrip_purp_mealsXaccess_time = mode_rail_l*trip_purp_meals*access_time,
    mode_rail_lXtrip_purp_medicalXaccess_time = mode_rail_l*trip_purp_medical*access_time,
    mode_rail_lXtrip_purp_otherXaccess_time = mode_rail_l*trip_purp_other*access_time,
    mode_rail_lXtrip_purp_schoolXaccess_time = mode_rail_l*trip_purp_school*access_time,
    mode_rail_lXtrip_purp_shoppingXaccess_time = mode_rail_l*trip_purp_shopping*access_time,
    mode_rail_lXtrip_purp_socialXaccess_time = mode_rail_l*trip_purp_social*access_time,
    mode_rail_lXtrip_purp_transp_someoneXaccess_time = mode_rail_l*trip_purp_transp_someone*access_time,
    mode_rail_lXtrip_purp_workXaccess_time = mode_rail_l*trip_purp_work*access_time,
    
    # trip purpose and mode specific wait  time
    mode_rail_lXtrip_purp_homeXwait_time = mode_rail_l*trip_purp_home*wait_time,
    mode_rail_lXtrip_purp_mealsXwait_time = mode_rail_l*trip_purp_meals*wait_time,
    mode_rail_lXtrip_purp_medicalXwait_time = mode_rail_l*trip_purp_medical*wait_time,
    mode_rail_lXtrip_purp_otherXwait_time = mode_rail_l*trip_purp_other*wait_time,
    mode_rail_lXtrip_purp_schoolXwait_time = mode_rail_l*trip_purp_school*wait_time,
    mode_rail_lXtrip_purp_shoppingXwait_time = mode_rail_l*trip_purp_shopping*wait_time,
    mode_rail_lXtrip_purp_socialXwait_time = mode_rail_l*trip_purp_social*wait_time,
    mode_rail_lXtrip_purp_transp_someoneXwait_time = mode_rail_l*trip_purp_transp_someone*wait_time,
    mode_rail_lXtrip_purp_workXwait_time = mode_rail_l*trip_purp_work*wait_time) %>%
  mutate( # generate transit-rail_c related variables
    # mode specific constant
    mode_rail_c = mode_rail_c,
    
    mode_rail_cXtrip_s = mode_rail_c*trip_dist_bin_s,
    mode_rail_cXtrip_m = mode_rail_c*trip_dist_bin_m,
    mode_rail_cXtrip_l = mode_rail_c*trip_dist_bin_l,
    mode_rail_cXtrip_ll = mode_rail_c*trip_dist_bin_ll,
    
    
    
    # mode specific travel time
    mode_rail_cXinv_time = mode_rail_c*inv_time,
    mode_rail_cXaccess_time = mode_rail_c*access_time,
    mode_rail_cXwait_time = mode_rail_c*wait_time,
    
    #trip purpose and mode specific constant
    mode_rail_cXtrip_purp_home = mode_rail_c*trip_purp_home,
    mode_rail_cXtrip_purp_meals = mode_rail_c*trip_purp_meals,
    mode_rail_cXtrip_purp_medical = mode_rail_c*trip_purp_medical,
    mode_rail_cXtrip_purp_other = mode_rail_c*trip_purp_other,
    mode_rail_cXtrip_purp_school = mode_rail_c*trip_purp_school,
    mode_rail_cXtrip_purp_shopping = mode_rail_c*trip_purp_shopping,
    mode_rail_cXtrip_purp_social = mode_rail_c*trip_purp_social,
    mode_rail_cXtrip_purp_transp_someone = mode_rail_c*trip_purp_transp_someone,
    mode_rail_cXtrip_purp_work = mode_rail_c*trip_purp_work,
    
    # trip purpose and mode specific inv travel time
    mode_rail_cXtrip_purp_homeXinv_time = mode_rail_c*trip_purp_home*inv_time,
    mode_rail_cXtrip_purp_mealsXinv_time = mode_rail_c*trip_purp_meals*inv_time,
    mode_rail_cXtrip_purp_medicalXinv_time = mode_rail_c*trip_purp_medical*inv_time,
    mode_rail_cXtrip_purp_otherXinv_time = mode_rail_c*trip_purp_other*inv_time,
    mode_rail_cXtrip_purp_schoolXinv_time = mode_rail_c*trip_purp_school*inv_time,
    mode_rail_cXtrip_purp_shoppingXinv_time = mode_rail_c*trip_purp_shopping*inv_time,
    mode_rail_cXtrip_purp_socialXinv_time = mode_rail_c*trip_purp_social*inv_time,
    mode_rail_cXtrip_purp_transp_someoneXinv_time = mode_rail_c*trip_purp_transp_someone*inv_time,
    mode_rail_cXtrip_purp_workXinv_time = mode_rail_c*trip_purp_work*inv_time,
    
    # trip purpose and mode specific access  time
    mode_rail_cXtrip_purp_homeXaccess_time = mode_rail_c*trip_purp_home*access_time,
    mode_rail_cXtrip_purp_mealsXaccess_time = mode_rail_c*trip_purp_meals*access_time,
    mode_rail_cXtrip_purp_medicalXaccess_time = mode_rail_c*trip_purp_medical*access_time,
    mode_rail_cXtrip_purp_otherXaccess_time = mode_rail_c*trip_purp_other*access_time,
    mode_rail_cXtrip_purp_schoolXaccess_time = mode_rail_c*trip_purp_school*access_time,
    mode_rail_cXtrip_purp_shoppingXaccess_time = mode_rail_c*trip_purp_shopping*access_time,
    mode_rail_cXtrip_purp_socialXaccess_time = mode_rail_c*trip_purp_social*access_time,
    mode_rail_cXtrip_purp_transp_someoneXaccess_time = mode_rail_c*trip_purp_transp_someone*access_time,
    mode_rail_cXtrip_purp_workXaccess_time = mode_rail_c*trip_purp_work*access_time,
    
    # trip purpose and mode specific wait  time
    mode_rail_cXtrip_purp_homeXwait_time = mode_rail_c*trip_purp_home*wait_time,
    mode_rail_cXtrip_purp_mealsXwait_time = mode_rail_c*trip_purp_meals*wait_time,
    mode_rail_cXtrip_purp_medicalXwait_time = mode_rail_c*trip_purp_medical*wait_time,
    mode_rail_cXtrip_purp_otherXwait_time = mode_rail_c*trip_purp_other*wait_time,
    mode_rail_cXtrip_purp_schoolXwait_time = mode_rail_c*trip_purp_school*wait_time,
    mode_rail_cXtrip_purp_shoppingXwait_time = mode_rail_c*trip_purp_shopping*wait_time,
    mode_rail_cXtrip_purp_socialXwait_time = mode_rail_c*trip_purp_social*wait_time,
    mode_rail_cXtrip_purp_transp_someoneXwait_time = mode_rail_c*trip_purp_transp_someone*wait_time,
    mode_rail_cXtrip_purp_workXwait_time = mode_rail_c*trip_purp_work*wait_time)

TravelMode.modmat = tmp
save(TravelMode.modmat, file = file.path(outdir, 'NHTS/TravelMode.modmat.National_transgeo.RData'))

# get a smaller set
tmp = TravelMode %>%
  rename(trip_purp = trip_purpose)%>%
  fastDummies::dummy_cols(select_columns = c('trip_purp','mode','trip_dist_bin')) %>%
  mutate(# generate covariates for auto nest modes (current only one)
    
    # mode specific constant
    mode_hv = mode_hv,
    
    # mode specific constant interacting with distance bins
    mode_hvXtrip_s = mode_hv*trip_dist_bin_s,
    mode_hvXtrip_m = mode_hv*trip_dist_bin_m,
    mode_hvXtrip_l = mode_hv*trip_dist_bin_l,
    mode_hvXtrip_ll = mode_hv*trip_dist_bin_ll,
    
    # mode specific travel time
    mode_hvXinv_time = mode_hv*inv_time
    
 ) %>%
  mutate( # generate walk mode variables
    # mode specific constant
    mode_walk = mode_walk,
    # mode specific constant interacting with distance bins
    mode_walkXtrip_s = mode_walk*trip_dist_bin_s,
    mode_walkXtrip_m = mode_walk*trip_dist_bin_m,
    mode_walkXtrip_l = mode_walk*trip_dist_bin_l,
    mode_walkXtrip_ll = mode_walk*trip_dist_bin_ll,
    
    # mode specific travel time
    mode_walkXinv_time = mode_walk*inv_time
    
  ) %>%
  mutate( # generate bike related variables
    # mode specific constant
    mode_bike = mode_bike,
    # mode specific constant interacting with distance bins
    mode_bikeXtrip_s = mode_bike*trip_dist_bin_s,
    mode_bikeXtrip_m = mode_bike*trip_dist_bin_m,
    mode_bikeXtrip_l = mode_bike*trip_dist_bin_l,
    mode_bikeXtrip_ll = mode_bike*trip_dist_bin_ll,
    
    
    # mode specific travel time
    mode_bikeXinv_time = mode_bike*inv_time
    
     ) %>%
  mutate( # generate transit-bus related variables
    # mode specific constant
    mode_bus = mode_bus,
    
    mode_busXtrip_s = mode_bus*trip_dist_bin_s,
    mode_busXtrip_m = mode_bus*trip_dist_bin_m,
    mode_busXtrip_l = mode_bus*trip_dist_bin_l,
    mode_busXtrip_ll = mode_bus*trip_dist_bin_ll,
    
    
    # mode specific travel time
    mode_busXinv_time = mode_bus*inv_time,
    mode_busXaccess_time = mode_bus*access_time,
    mode_busXwait_time = mode_bus*wait_time) %>%
  mutate( # generate transit-taxi related variables
    # mode specific constant
    mode_taxi = mode_taxi,

    mode_taxiXtrip_s = mode_taxi*trip_dist_bin_s,
    mode_taxiXtrip_m = mode_taxi*trip_dist_bin_m,
    mode_taxiXtrip_l = mode_taxi*trip_dist_bin_l,
    mode_taxiXtrip_ll = mode_taxi*trip_dist_bin_ll,
    
    
    # mode specific travel time
    mode_taxiXinv_time = mode_taxi*inv_time,
    mode_taxiXaccess_time = mode_taxi*access_time,
    mode_taxiXwait_time = mode_taxi*wait_time) %>%
  mutate( # generate transit-rail_l related variables
    # mode specific constant
    mode_rail_l = mode_rail_l,
 
    mode_rail_lXtrip_s = mode_rail_l*trip_dist_bin_s,
    mode_rail_lXtrip_m = mode_rail_l*trip_dist_bin_m,
    mode_rail_lXtrip_l = mode_rail_l*trip_dist_bin_l,
    mode_rail_lXtrip_ll = mode_rail_l*trip_dist_bin_ll,
    
    # mode specific travel time
    mode_rail_lXinv_time = mode_rail_l*inv_time,
    mode_rail_lXaccess_time = mode_rail_l*access_time,
    mode_rail_lXwait_time = mode_rail_l*wait_time) %>%
  mutate( # generate transit-rail_c related variables
    # mode specific constant
    mode_rail_c = mode_rail_c,

    mode_rail_cXtrip_s = mode_rail_c*trip_dist_bin_s,
    mode_rail_cXtrip_m = mode_rail_c*trip_dist_bin_m,
    mode_rail_cXtrip_l = mode_rail_c*trip_dist_bin_l,
    mode_rail_cXtrip_ll = mode_rail_c*trip_dist_bin_ll,
    
    
    # mode specific travel time
    mode_rail_cXinv_time = mode_rail_c*inv_time,
    mode_rail_cXaccess_time = mode_rail_c*access_time,
    mode_rail_cXwait_time = mode_rail_c*wait_time)

TravelMode.modmat = tmp
save(TravelMode.modmat, file = file.path(outdir, 'NHTS/TravelMode.modmat.National.small_transgeo.RData'))

load(file = file.path(outdir, 'NHTS/TravelMode.modmat.National.small_transgeo.RData'))







