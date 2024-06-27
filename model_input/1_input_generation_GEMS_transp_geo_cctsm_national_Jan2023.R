# Generating model inputs for GEMS
# NATALIE POPOVICH
# BERKELEY LAB
# SEP 10 2020
# UPDATED MARCH 29 2021: include TransitionMatrix, Microtype Assignment, and  AvgTripLengths
# Updated May 11 2021: transition matrix using ROUTED NHTS trips and correct tract sequencing 
# Updated May 14 2021: transition matrix using ROUTED NHTS trips and correct tract sequencing at CCST level
# Updated May 19 2021: changing trip generation rates to hourly time bins
# Updated Aug 26 2021: adding pooled betas to the same file as the non-pooled betas for the mode choice model
# Updated Aug 28 2021: Michelle
# Updated Sep 13 2021: Natalie added CCST level data with imputed clusters
# Udpdated Nov 5 2021: Natalie added some of the Task 2 output files, restructured variables and file names
# Updated Nov 11 2021: changing file structure a bit and including additional mode availability layer, combining system costs into 
    #mode-specific cost files
# Updated Dec 8 2021: adding bikeshare system density to bike.csv file
# Updated Dec 9 2021: adding Microtype Mode Availability layer
# Updated Dec 14 2021: added average distance between rail stations, updated bikeshare costs to be total per bike per year, removed extra variables 
              #from the bus and rail files
# Updated Jan 11 2022: using mode speeds from the NHTS for rail, walk, and bike
# Updated March 10 2022: fixed transition matrix so that we don't have A_1 to A_1 transition etc... (diagonals should be zero)
# changed microtype assignment to accomodate new transition matrix approach
# changed avg thru lengths  to accomodate new transition matrix approach

# Xiaodan Updates (Jan 2023):
# 1. update distance bin 2. add geotypes to demand generation (trip rate) 
# 3. update bike share 4. commented codes that are no longer needed

# Xiaodan Updates (Mar 2023):
# 1. adding job count into activity density. 
# 2. removing crash and ghg costs from externality
######################

#files exported:
#AvgTripLengths.csv  --> no longer needed
#DistanceBins.csv  --> 8 bins
#DistanceDistribution.csv --> 8 bins
#FreightDemand.csv
#MicrotypeAssignment.csv --> no longer needed
#Microtypes.csv --> no longer needed
#ModeAvailability.csv
#OriginDestination.csv 
#Population.csv 
#PopulationGroups.csv
#TimePeriods.csv 
#TransitionMatrix.csv
#TripGeneration.csv --> add geotype
#MicrotypeModeAvailability.csv
#TripPurposes.csv 

#MODES
#walk.csv
#auto.csv
#rail.csv
#bus.csv
#bike.csv
#ridehail.csv


# Set working directory and subdirectories
rm(list=ls())
mywd = "C:/FHWA/For FHWA folks/trip_generation_and_gems_inputs"
#setwd() # 
#setwd("D:/Box Sync/FHWA/Task3/") # Michelle
setwd(mywd) # Xiaodan

# datadir <- "./data/Raw"
cleandir <- "./GEMS_inputs"
resultsdir <- './results/TransportGeography'
modeldir <-  './results/TransportGeography/ModelInputs'
modedir <- './results/TransportGeography/ModelInputs/modes'

# library(sqldf) # use sql commands for csv files 
library(tidyverse)
library(data.table)
library(stats)
library(openxlsx)

########
# LOAD DATASETS
#########

# mode speeds from NHTS
speeds = read.xlsx(file.path(cleandir, "mode_speeds.xlsx"))

# NHTS TRIPS
#trips <- fread(file.path(datadir, "nhts_no_ids.csv")) #cleaned from the NHTS 
# generated from 1_nhts_data_anonymize.R on the server
trips = fread(file.path(cleandir, "nhts_no_ids_1hrtimebins_with_imputation.csv")) %>% # smaller time bins
  rename(PERSONID = PERSON_ID)

#summarize the missing values
trips %>% # Here none of the microtypeIDs are missing, but some of the Geotypes are blank
  select(everything()) %>%
  summarise_all(funs(sum(is.na(.))))

trips[trips==""]<-NA #825 geotypes are still missing for the ODs, but showed up as blanks

# impute origin and destination geotypes using home geotype
trips = trips %>% 
  dplyr::mutate(o_geotype = case_when(
                    !is.na(o_geotype) ~ o_geotype, 
                    TRUE ~ h_geotype),
          d_geotype = case_when(
                      !is.na(d_geotype) ~ d_geotype, 
                      TRUE ~ h_geotype))

# NHTS POPULATION ESTIMATES 
# Updated Nov 11 2021: new population classes: income, vehicle ownership, senior
pop = fread(file.path(cleandir, "nhts_user_classes_inc_veh_sr.csv")) %>% 
  select(PERSONID, HOUSEID, PopulationGroupID, WTPERFIN, person_indx) %>%
  mutate(PopulationUS = as.integer(sum(WTPERFIN)))

# Transport Geo cluster labels ( this might be called "type" below, should change to "labels")
# NOTE SEP 10 2021: this now uses updated CCSTS using imputed micro-geotypes
labels = read.csv(file.path(cleandir, "ccst_geoid_key_transp_geo_with_imputation.csv")) %>%
  mutate(GEOID = str_pad(as.character(GEOID), width = 11, side = "left", pad = "0")) %>%
  unite("MicrotypeID", geotype, microtype, sep = "_", remove = F) %>%
  select(spatial_id, FID, GEOID, MicrotypeID, microtype, geotype, cbsa)

tract_labels = read.xlsx(file.path(cleandir, "micro_geo_results_with_imputation_volpe.xlsx")) %>%
  mutate(GEOID = str_pad(as.character(GEOID), width = 11, side = "left", pad = "0"))

# Output from mode choice model by user groups
# Mode costs (operating and user)
modecost = read.csv(file.path(cleandir, "modecost.csv"))  %>%
  left_join(tract_labels %>% select(MicrotypeID, geotype)) %>%
  unique()

#mode choice coefficients
popgroups = read.csv(file.path(cleandir, "populationgroups.csv"))

# Opportunities per tract (we're not really using this right now)
ops <- fread(file.path(cleandir, "opportunity_counts_tract.csv")) %>%
  mutate(GEOID = str_pad(as.character(GEOID), width = 11, side = "left", pad = "0"))

ops_jobs <- fread(file.path(cleandir, "wac_tract_2017.csv")) %>%
  mutate(GEOID = str_pad(as.character(trct), width = 11, side = "left", pad = "0"))

# Note: updated NHTS data May 14 2021 at CCST level
# NOTE SEP 13 2021: updated to included imputed clusters
lengths = fread(file.path(cleandir, "nhts_thru_lengths_ordered_ccst_transp_geo_with_imputation.csv")) %>%
  mutate(thru_geoid = str_pad(as.character(thru_geoid), width = 11, side = "left", pad = "0"),
          o_geoid = str_pad(as.character(o_geoid), width = 11, side = "left", pad = "0"),
          d_geoid = str_pad(as.character(d_geoid), width = 11, side = "left", pad = "0"))

# Read in NHTS OD pairs and person weights
# NOTE Sep 13 2021: this file is okay becuase it doens't have geotype info in it
nhts= fread(file.path(cleandir, "nhts_od_pairs_2017.csv")) %>%
  mutate(o_geoid = str_pad(as.character(o_geoid), width = 11, side = "left", pad = "0"),
        d_geoid = str_pad(as.character(d_geoid), width = 11, side = "left", pad = "0"))

#rail station distance
raildistance = read.csv(file.path(cleandir, "raildistance.csv"))

# #import external costs to generate mode costs
ext <- fread(file.path(cleandir, "external_costs_mode_microtype_transp_geo.csv")) %>%
  mutate(PerMileCost = cost_air+ cost_noise) %>%
  select(MicrotypeID, PerMileCost, Mode) 

# Road network lenght for each tract 
row = read.csv(file.path(cleandir, "row.csv")) %>%
     mutate(GEOID = str_pad(as.character(tract), width = 11, side = "left", pad = "0"))
             
row2 = fread(file.path(cleandir, "row_no_sample.csv")) %>%
  mutate(GEOID = str_pad(as.character(tract), width = 11, side = "left", pad = "0"))

#Transit service with all modes binned
service = read.csv(file.path(cleandir, "TransitService.csv")) %>% # note these are the transport geography results from Task 2 folder
  rename(VehicleRevenueHours = veh_revenue_hours,
         VehicleRevenueMiles = veh_revenue_miles,
         DirectionalRouteMiles= dir_route_miles, 
         Mode = mode_group) 

# Transit layer
transitlayer = read.csv(file.path(cleandir, "modeaccessibility.csv")) %>%
  mutate(GEOID = str_pad(as.character(geoid), width = 11, side = "left", pad = "0")) %>%
    select(-X, -geoid)

bikeshare = read.csv(file.path(cleandir, "bikeshare_national.csv"))%>%
  select(-Mode)

bikeshare$StationDensity[is.na(bikeshare$StationDensity)] <- 0

# Read in NHTS OD pairs and person weights and CCST level
# IGNORE THIS FILE, NOT SURE WHEN IT WAS USED
#nhts= fread(file.path(datadir, "nhts_od_pairs_2017_ccst_transp_geo.csv")) 

#road network costs
network = read.csv(file.path(cleandir, "RoadNetworkCosts.csv"))

###############################
# DATA CLEANING
##############################
#summarize the missing values
trips %>%
  select(everything()) %>%
  summarise_all(funs(sum(is.na(.)))) # nothing missing here except median income (which we don't use in this file)

# Making naming convention match Zach's input files
trips = trips %>% 
  #select(-spatial_id) %>%
  unite('OriginMicrotypeID', c(o_geotype,o_microtype), sep = "_", remove = T) %>%
  unite('HomeMicrotypeID', c(h_geotype,h_microtype), sep = "_", remove = T) %>%
  unite('DestinationMicrotypeID', c(d_geotype,d_microtype), sep = "_", remove = T) %>%
  # FOCUSING ONLY ON 6 TRIP PURPOSES: home, work, school, medical, social (meals/rec) + other
  mutate(TripPurposeID = case_when(      
    trip_purpose %in% c("meals", "social", "shopping") ~ 'leisure',
    trip_purpose %in% c("other", "transp_someone") ~ 'other',
    TRUE ~ trip_purpose)) %>% 
  rename(TimePeriodID = start_time_bin) #,
#population_tract = population) # I think we want to ignore this since we can't link HHIDs to GEOIDs in this file

########################
# DistanceBins.csv
# variables: DistanceBinID	MeanDistanceInMiles
#########################
# trips = trips %>%
#   filter(!is.na(trpmiles))%>% #there should be no missing values anyway, but leave this just in case
#   mutate(DistanceBinID= case_when(
#     trpmiles <=1.3  ~ 'short',
#     trpmiles >1.3 & trpmiles<=3 ~ 'medium',
#     trpmiles>3 & trpmiles <=8 ~ 'long',
#     TRUE ~ 'xlong'))

trips = trips %>%
  filter(!is.na(trpmiles))%>% #there should be no missing values anyway, but leave this just in case
  mutate(DistanceBinID= case_when(
    trpmiles <=1  ~ 'bin1',
    trpmiles >1 & trpmiles<=2 ~ 'bin2',
    trpmiles >2 & trpmiles<=4 ~ 'bin3',
    trpmiles >4 & trpmiles<=8 ~ 'bin4',
    trpmiles >8 & trpmiles<=15 ~ 'bin5',
    trpmiles >15 & trpmiles<=20 ~ 'bin6',
    trpmiles>20 & trpmiles <=35 ~ 'bin7',
    TRUE ~ 'bin8'))

#NOTE: Updated Sep 14, 2021 to include person weights, though it appears not to have changed anything
# trips <- trips %>% filter(trpmiles <= quantile(trpmiles, 0.975))

distbin <- trips %>% 
  select(OriginMicrotypeID, DistanceBinID, trpmiles, wtperfin) %>% 
  group_by(OriginMicrotypeID, DistanceBinID) %>%
  mutate(MeanDistanceInMiles = weighted.mean(trpmiles, weight = wtperfin)) %>%
  select(-trpmiles, - wtperfin) %>% 
  unique()

fwrite(distbin, file.path(modeldir, "DistanceBins.csv"), row.names = F)

############################
# Population.csv
# variables: MicrotypeID, PopulationGroupID, Population 
##############################
#but let's check that merge worked correctly
test = pop %>% select(PERSONID, HOUSEID, person_indx) %>% distinct()
test2 = trips %>% select(PERSONID, HOUSEID) %>% distinct() 

#NOTE: trip weights are equal to 365 times person weights to generate one year of travel
trips = trips %>%
  right_join(pop, by = c("PERSONID", "HOUSEID")) # this file makes sense now because a number of people don't travel

trips %>% select(everything()) %>% summarise_all(funs(sum(is.na(.)))) #check missing geotypes again here
#no longer any trips without people files, only people who don't take trips
# NOTE: not missing any geotypes yet here

# here is where we lose a lot of trips because many people do not travel. 
# for accurate populatino estimates, I think we need to assign home geotypes to the person file first 

#generate total population from the person weights by population Group and microtype
pop <- trips %>% 
  select(HomeMicrotypeID, PopulationGroupID, WTPERFIN, person_indx) %>% 
  distinct() %>%
  group_by(HomeMicrotypeID, PopulationGroupID) %>%
  mutate(Population = as.integer(sum(WTPERFIN))) %>%
  rename(MicrotypeID = HomeMicrotypeID) %>%
  select(MicrotypeID, Population, PopulationGroupID) %>% 
  distinct() %>%
  ungroup()

fwrite(pop, file.path(modeldir, "Population.csv"), row.names = F) 

pop <- pop %>%
  separate(MicrotypeID, c("HomeGeotypeID", "HomeMicrotypeID"), remove = FALSE)

pop <- pop %>% # note this population only includes people who traveled in the NHTS
  group_by(PopulationGroupID) %>%
  mutate(PopulationPopGroup = sum(Population, na.rm = T)) %>% # generate total population of the group regardless of microtype
  ungroup() %>%
  group_by(HomeGeotypeID, PopulationGroupID) %>%
  mutate(PopulationPopGroupGT = sum(Population, na.rm = T)) %>% # generate total population of the group regardless of microtype
  ungroup() %>%
  group_by(MicrotypeID) %>%
  mutate(PopulationMicrotypeID = sum(Population, na.rm = T))

########################
# TimePeriods.csv
# variables: TimePeriodID	DurationInHours
#########################
time <- trips %>% 
  ungroup() %>%
  select(TimePeriodID) %>% 
  distinct() %>%
  mutate(DurationInHours = 1) 

fwrite(time, file.path(modeldir,"TimePeriods.csv"), row.names = F)

################
# TripPurposes.csv
# variables: TripPurposeID
################
purpose <- trips %>% 
  ungroup() %>%
  select(TripPurposeID) %>% 
  distinct() %>%
  filter(!is.na(TripPurposeID))

fwrite(purpose, file.path(modeldir, "TripPurposes.csv"), row.names = F)

################
# Microtypes.csv
# variables: MicrotypeID
################
# type <- pop %>% 
#   ungroup() %>%
#   select(MicrotypeID) %>% 
#   distinct() %>%
#   arrange(MicrotypeID) %>%
#   filter(!is.na(MicrotypeID))
# 
# fwrite(type, file.path(modeldir, "Microtypes.csv"), row.names = F)

#############
# Generating total number of trips so that they can be apportioned by distance bins and trip rates
# TripGeneration.csv
# variables: TimePeriodID	PopulationGroupID	TripPurposeID	TripGenerationRatePerHour
#########################

# the person weights should correspond to each trip purpose and time bin 
# note that these weights only include people who took trips!!
trips <- trips %>%
  separate(HomeMicrotypeID, c("HomeGeotype", "HomeMicrotype"), remove = FALSE)

rates <- trips %>%
  filter(!is.na(HomeMicrotypeID)) %>%
  group_by(HomeGeotype, TimePeriodID,PopulationGroupID,TripPurposeID) %>%
  add_tally(name = "total_trips_purp_time", wt = WTPERFIN) %>% 
  merge(time, by = "TimePeriodID") %>%
  mutate(MicrotypeID = HomeMicrotypeID) %>%  # these rates are independent of microtype, and just based on population groups
  #bring population of each population group back in 
  # these need to be adjusted to include only the population of people who took trips
  left_join(pop, by = c("MicrotypeID", "PopulationGroupID")) %>% # Nov 17 21: check if this should be left or right join
  # divide by population of home micro-geo pair and number of hours in each time bin
  mutate(TripGenerationRatePerHour = total_trips_purp_time/(PopulationPopGroupGT*DurationInHours)) %>%
  select(HomeGeotype, TimePeriodID, PopulationGroupID, TripPurposeID, TripGenerationRatePerHour) %>%
  distinct() %>%
  filter(!is.na(TripPurposeID))

head(rates)
fwrite(rates, file.path(modeldir, "TripGeneration.csv"), row.names = F)

########################
# PopulationGroups.csv
# variables: PopulationGroupID	TripPurposeID Mode Intercept 
# BetaWaitAccessTime BetaTravelTime BikeShare_Bike BetaMonetaryCost 
# BetaTravelTime_Pooled  BetaWaitAccessTime_Pooled BetaMonetaryCost_Pooled Intercept_Pooled BikeShare_Bike_Pooled

# uses output from mode choice model
#########################
# add seniors to user classes
seniors = popgroups %>%
  mutate(PopulationGroupID = case_when(
    PopulationGroupID == "HighIncNoVeh" ~ "HighIncNoVehSenior",
    PopulationGroupID == "HighIncVeh" ~ "HighIncVehSenior",
    PopulationGroupID == "LowIncNoVeh" ~ "LowIncNoVehSenior", 
    PopulationGroupID == "LowIncVeh" ~ "LowIncVehSenior"))

popgroups = rbind(popgroups, seniors)

# add trip purposes, even though for now the don't vary 
popgroups = popgroups %>%
  mutate(TripPurposeID = "work")

school = popgroups %>% 
  mutate(TripPurposeID = "school")

home = popgroups %>% 
  mutate(TripPurposeID = "home")

leisure = popgroups %>%
  mutate(TripPurposeID = "leisure")

other = popgroups %>% 
  mutate(TripPurposeID = "other")

medical= popgroups %>% 
  mutate(TripPurposeID = "medical")

popgroups = rbind(popgroups, home) %>% 
  rbind(leisure) %>% 
  rbind(other) %>% 
  rbind(medical) %>% 
  unique()

fwrite(popgroups, file.path(modeldir, "PopulationGroups.csv"), row.names = F)

########################
# OriginDestination.csv
# variables: HomeMicrotypeID	TimePeriodID	PopulationGroupID	 TripPurposeID	OriginMicrotypeID	DestinationMicrotypeID	Portion
# Portion is distributed across destination microtypes, all else equal
#########################
ods = trips %>% 
  select(HomeMicrotypeID, TimePeriodID,	PopulationGroupID, TripPurposeID,	OriginMicrotypeID,	DestinationMicrotypeID,  wtperfin) %>%
  group_by(HomeMicrotypeID, TimePeriodID,	PopulationGroupID, TripPurposeID,	OriginMicrotypeID,	DestinationMicrotypeID) %>%
  add_tally(name = "total_trips_detailed", wt = wtperfin) %>% # person weights
  group_by(HomeMicrotypeID, TimePeriodID,	PopulationGroupID, TripPurposeID) %>% # assign proportion of trips by OD pairs
  add_tally(name = "denominator", wt = wtperfin) %>% # person weights
  mutate(Portion = total_trips_detailed/denominator) %>% 
  select(HomeMicrotypeID, TimePeriodID,	PopulationGroupID, TripPurposeID,	OriginMicrotypeID,DestinationMicrotypeID, Portion) %>%
  distinct() %>%
  filter(!is.na(OriginMicrotypeID))

fwrite(ods, file.path(modeldir, "OriginDestination.csv"), row.names = F)

########################
# DistanceDistribution.csv
# variables: TripPurposeID	OriginMicrotypeID	DestinationMicrotypeID	DistanceBinID	Portion
# assigns total number of trips for the same OD pairs to distance bins
#########################
# generate total number of trips between each OD pair for each home microtyope 
dd = trips %>% 
  group_by(OriginMicrotypeID, DestinationMicrotypeID, TripPurposeID, DistanceBinID) %>% 
  add_tally(name = "trips_bin", wt = wtperfin) %>% # person weights 
  group_by(OriginMicrotypeID, DestinationMicrotypeID, TripPurposeID) %>% 
  add_tally(name = "trips_all", wt = wtperfin) %>%
  mutate(Portion = trips_bin/trips_all) %>%
  ungroup() %>%
  select(TripPurposeID, OriginMicrotypeID, DestinationMicrotypeID, DistanceBinID, Portion) %>%
  distinct() %>%
  arrange(OriginMicrotypeID, DestinationMicrotypeID, TripPurposeID)

fwrite(dd, file.path(modeldir, "DistanceDistribution.csv"), row.names = F)

########################
# MicrotypeModeAvailability.csv
# Describes the percent of each microtype geotype that has each transit layer
# variables: MicrotypeID TransitLayer Portion
# Portion is distributed across types of transit availability 
#########################
layers = transitlayer %>%
  left_join(tract_labels %>% select(GEOID, MicrotypeID)) %>% 
  mutate(TransitLayer = case_when( 
           bike ==1 & bus == 1 & rail == 1 ~"bus_rail_bike",
           bike ==0 & bus == 1 & rail == 1 ~ "bus_rail_nobike",
           bike ==1 & bus == 1 & rail == 0 ~ "bus_norail_bike", 
           bike == 0 & bus == 1 & rail == 0 ~ "bus_norail_nobike",
           bike ==1 & bus == 0 & rail == 1 ~ "nobus_rail_bike",
           bike ==1 & bus == 0 & rail == 0 ~ "nobus_norail_bike",
           bike ==0 & bus == 0 & rail == 1 ~ "nobus_rail_nobike",
           bike ==0 & bus == 0 & rail == 0 ~ "nobus_norail_nobike"))  %>%
  group_by(MicrotypeID) %>% # calculate percent of transitlayers in each micro-geotype
  add_tally() %>% 
  group_by(MicrotypeID, TransitLayer) %>%
  add_tally() %>% 
  ungroup() %>%
  mutate(Portion = nn/n) %>%
  select(MicrotypeID, TransitLayer, Portion) %>%
  unique() %>%
  arrange(MicrotypeID)
  
write.csv(layers, file.path(modeldir, "MicrotypeModeAvailability.csv"), row.names = F)

########################
# ModeAvailability.csv
# variables: OriginMicrotypeID	DestinationMicrotypeID TransitLayer Portion
# Portion is distributed across types of transit availability 
# Apportions total number of trips from each OD microtype pair across transit availability layers
# Note this is at the TRACT level not CCST level since transit layers may vary within a CCST 
#########################
# create data tables for smaller file sizes
nhts =as.data.table(nhts)  # each row is one trip ID

avail = as.data.table(lengths) %>% #distance traveled by each trip in each tract
  filter(!is.na(o_geoid)) %>%
  rename(GEOID = thru_geoid) %>%
  select(-o_fid, -d_fid, -thru_length_ccst, -thru_fid)

avail = nhts %>%
  left_join(avail) %>%
  unique() %>%
  group_by(trip_indx) %>%
  mutate(trip_length_miles = distance/1609.344) %>%
  ungroup() %>%
  select(trip_indx, trip_length_miles, o_geoid, d_geoid, wtperfin) %>% 
  distinct()
head(avail)

# Now need to calculate portion of each trip that went through each level of transit availability
# in this case, we don't care about the "thru" geoids or distances, only the OD and total trip legnth
# need to merge origin and destination microtypes to each o_geoid and d_geoid
avail2 = avail %>%
  rename(GEOID = o_geoid) %>% 
  left_join(transitlayer) %>% # merge origin tract with transit layer
  rename(o_bike = bike, o_bus =bus, o_rail = rail) %>%
  left_join(tract_labels %>% select(GEOID, MicrotypeID)) %>% # merge origin tract with correct label
  rename(OriginMicrotypeID = MicrotypeID) %>%
  select(-GEOID) %>%
  rename(GEOID = d_geoid) %>% 
  left_join(transitlayer) %>% # merge origin tract with transit layer
  rename(d_bike = bike, d_bus =bus, d_rail = rail) %>%
  # merge origin tract with transit layer
  left_join(tract_labels %>% select(GEOID, MicrotypeID)) %>% # merge dest tract with correct label
  rename(DestinationMicrotypeID = MicrotypeID) 

head(avail2)

#assign transit availability based on least available mode
avail2 = avail2 %>%
  mutate(bike = case_when(o_bike ==1 & d_bike == 1 ~ 1,    TRUE ~ 0),
         bus = case_when(o_bus ==1 & d_bus == 1 ~ 1,    TRUE ~ 0),
         rail = case_when(o_rail ==1 & d_rail == 1 ~ 1,    TRUE ~ 0),
         TransitLayer = case_when( # generate transit layer based on combined OD transit availability
           bike ==1 & bus == 1 & rail == 1 ~"bus_rail_bike",
           bike ==0 & bus == 1 & rail == 1 ~ "bus_rail_nobike",
           bike ==1 & bus == 1 & rail == 0 ~ "bus_norail_bike", 
           bike == 0 & bus == 1 & rail == 0 ~ "bus_norail_nobike",
           bike ==1 & bus == 0 & rail == 1 ~ "nobus_rail_bike",
           bike ==1 & bus == 0 & rail == 0 ~ "nobus_norail_bike",
           bike ==0 & bus == 0 & rail == 1 ~ "nobus_rail_nobike",
           bike ==0 & bus == 0 & rail == 0 ~ "nobus_norail_nobike"))

head(avail2)
# assign portion of total VMT by origin destination microtype
avail2 = avail2 %>%
  group_by(OriginMicrotypeID, DestinationMicrotypeID) %>%
  mutate(DEN = sum(trip_length_miles*wtperfin, na.rm = T)) %>%
  group_by(OriginMicrotypeID, DestinationMicrotypeID, TransitLayer) %>%
  mutate(Portion = sum(wtperfin*trip_length_miles, na.rm = T)/DEN) %>%
  select(OriginMicrotypeID, DestinationMicrotypeID, TransitLayer, Portion) %>%
  unique()

fwrite(avail2, file.path(modeldir, "ModeAvailability.csv"), row.names = F)

###############
# MODE ATRRIBUTES 
##############

# convert speeds from miles per hour to meters per second
speeds = speeds %>%
  mutate(SpeedInMetersPerSecond = mean.speed.m_h/2.237) %>%
  rename(Mode = mode_chosen) %>%
  filter(start_time_bin == "morning_rush") %>%
  select(Mode, SpeedInMetersPerSecond) %>%
  unique()

#Note: these data are estimated somewhere in excel and so for now are manually constructed from the results of the 
# Task2 system costs calibration exercsie
# They also combined operating and capital costs since the regression is for total costs over the year on fleet size

# bike.csv
# Variables: MicrotypeID, PerStartCost, PerMinuteCost, DailyCapCostPerBike, DailyOpCostPerBike
# estimates from CapitalBikeShare system
##############
#NOTE: dec 14 2021 need to update to have only one variable cost per bike
bike = data.frame("geotype" = c("A", "B", "C", "D", "E" ,"F"), 
                  "CapCostPerBike" = 9000 , "DailyOpCostPerBike" = 3400/365,
                  "CapCostPerDock" = 4500 , "DailyOpCostPerDock" = 1700/365, 
                  "DocksPerBike" = 2, "DocksPerStation" = 15,
                  "PerStartCost" = 2.5, "PerMinuteCost" = .05,
                  "Mode" = "bike") %>% 
  right_join(tract_labels %>% select(MicrotypeID, geotype)) %>%
  unique() %>%
  select(MicrotypeID, everything(.), -geotype)  %>%
  left_join(bikeshare)  %>%
  left_join(speeds) %>%
  select(-Mode)

fwrite(bike, file.path(modedir, "bike.csv"), row.names = F)

################
# auto.csv
# variables: MicrotypeID	PerStartCost	PerEndCost	PerMileCost	VehicleSize
################
auto = modecost %>% 
  filter(Mode == "auto") %>%
  mutate(VehicleSize = NA) %>%
  select(MicrotypeID,	PerStartCost,	PerEndCost,	PerMileCost,	VehicleSize)

fwrite(auto, file.path(modedir, "auto.csv"), row.names = F)

########################
# bus.csv
# variables: MicrotypeID	DailyCostPerVehicle Headway	PerStartCost	SeniorFareDiscount VehicleCapacity	VehicleSize	StopSpacing
#########################
bus = modecost %>% 
  filter(Mode == "bus") %>%
  left_join(service %>% filter(Mode == "bus") ) %>%
  mutate(DailyCostPerVehicle = case_when(
    geotype == "A" ~2004,
    geotype == "B" ~ 1761, 
    geotype == "C" ~ 1370,
    geotype == "D" ~ 1292,
    geotype == "E" ~ 928,
    TRUE ~ 1529),
    SeniorFareDiscount = .5,
    Headway = NA,
    VehicleCapacity = NA,
    VehicleSize = NA,
    StopSpacingKM = NA) %>%
  select(-Mode, - PerMileCost, -PerMinuteCost, - PerEndCost) %>% 
  # need to distribute network length across microtypes within a geotype. 
  # for now, use the proportion of the population in each microgeo ID to weight
  left_join(pop %>% select(MicrotypeID, PopulationMicrotypeID)) %>%
  unique() %>%
  group_by(geotype) %>%
  mutate(weight = PopulationMicrotypeID/sum(PopulationMicrotypeID),
         VehicleRevenueHoursPerDay = VehicleRevenueHours*weight/365,
         VehicleRevenueMilesPerDay = VehicleRevenueMiles*weight/365, 
         DirectionalRouteMiles = DirectionalRouteMiles*weight) %>%
  ungroup() %>%
  select(-weight, -VehicleRevenueHours, -VehicleRevenueMiles) 

fwrite(bus, file.path(modedir,"bus.csv"), row.names = F)


########################
# rail.csv
# variables: MicrotypeID DailyCostPerVehicle	Headway	PerStartCost	PerMileCost	StopSpacing	SpeedInMetersPerSecond
#########################

#for now aggregate light rail and commuter rail into one rail mode for route lengths
railservice = service %>%
  filter(Mode != "bus") %>%
  group_by(geotype) %>%
  mutate(VehicleRevenueHours = sum(VehicleRevenueHours, na.rm = T),
         VehicleRevenueMiles = sum(VehicleRevenueMiles, na.rm = T),
         DirectionalRouteMiles = sum(DirectionalRouteMiles, na.rm = T)) %>%
  mutate(Mode = "rail") %>%
  unique()

# for now rail costs are estimate using heavy rail
rail = modecost %>% 
  filter(Mode == "rail") %>%
  left_join(railservice) %>%
  left_join(pop %>% select(MicrotypeID, PopulationMicrotypeID)) %>%
  left_join(speeds) %>%
  left_join(raildistance) %>%
  mutate(DailyCostPerVehicle = case_when(
    geotype == "A" ~ 3990,
    geotype == "B" ~ 3580, 
    TRUE ~ 4212),
    Headway = NA,
    VehicleCapacity = NA,
    VehicleSize = NA)  %>%
  rename(StopSpacingKM = Distance_km) %>%
  unique() %>%
  select(- Mode, - PerMileCost, -PerMinuteCost, - PerEndCost) %>%  
  group_by(geotype) %>%
  mutate(weight = PopulationMicrotypeID/sum(PopulationMicrotypeID),
         VehicleRevenueHoursPerDay = VehicleRevenueHours*weight/365,
         VehicleRevenueMilesPerDay = VehicleRevenueMiles*weight/365,
         DirectionalRouteMiles = DirectionalRouteMiles*weight) %>%
  ungroup() %>%
  select(- VehicleRevenueHours, -VehicleRevenueMiles, -weight, - geotype)

fwrite(rail, file.path(modedir, "rail.csv"), row.names = F)


########################
# walk.csv
# variables: MicrotypeID	SpeedInMetersPerSecond
#########################
walk = data.frame("geotype" = c("A", "B", "C", "D", "E" ,"F"),
                  Mode = "walk") %>%
  left_join(tract_labels %>% select(MicrotypeID, geotype)) %>%
  left_join(speeds) %>%
  unique() %>%
  select(MicrotypeID, SpeedInMetersPerSecond)

fwrite(walk, file.path(modedir, "walk.csv"), row.names = F)


########################
# ridehail.csv
# variables: MicrotypeID	PerStartCost	PerMileCost	VehicleCapacity	StopSpacing	SpeedInMetersPerSecond
#########################
ridehail = modecost %>%
  filter(Mode == "ridehail") %>%
  mutate(VehicleSize = NA,
         VehicleCapacity = NA) %>%
  select(MicrotypeID, PerStartCost, PerMileCost, PerMinuteCost, VehicleSize,	VehicleCapacity)

fwrite(ridehail, file.path(modedir, "ridehail.csv"), row.names = F)

########################
# EXTERNAL COSTS BY MODE 
# Externalities.csv
# variables: Mode MicrotypeID PerMileExtCost
# #############################
ext = ext %>% 
  rename(PerMileExtCost = PerMileCost)

fwrite(ext, file.path(modeldir, "Externalities.csv"), row.names = F)

# ########################
# # SubNetworks.csv
# # variables: SubnetworkID	MicrotypeID	ModesAllowed	LengthNetworkLaneMiles	vMax	Type
# #########################
subnet <- row %>%  # start with just the total length of ROW in each microtype 
  select(GEOID, lm_all_tract) %>% 
  distinct() %>%
  right_join(tract_labels, by = "GEOID") %>% 
  group_by(MicrotypeID) %>%
  mutate(LengthNetworkLaneMiles = sum(lm_all_tract, na.rm = T)) %>%
  select(MicrotypeID, LengthNetworkLaneMiles) %>%
  distinct() %>%
  arrange(MicrotypeID) %>%
  mutate(ModesAllowed = NA, 
         vMax = NA, 
         Type = NA)

fwrite(subnet, file.path(modeldir, "SubNetworks.csv"), row.names = F)

##########################
# RoadNetworkCosts.csv
# Note: currently generated in Task 2, but copied here for completeness
# Variables: Mode MicrotypeID, LaneDedicationPerLaneMile, ROWConstructionPerLaneMile
######################
write.csv(network, file.path(modeldir,"RoadNetworkCosts.csv"), row.names = F)

###################
# Accessibility Measure
# Need to calculate activity density per microtype
# MicrotypeID ActivityID OppsPerSqMi TripPurposeID
###################

ops_jobs = ops_jobs %>% 
  rename(num_jobs = jobs_total)

ops <- ops %>%
  merge(ops_jobs, by = "GEOID") %>%
  merge(tract_labels, by = "GEOID") %>%
  merge(row, by = "GEOID")  %>%
  select(contains("num"), aland, MicrotypeID, -num_jrcollege) %>%
  group_by(MicrotypeID) %>%
  mutate_at(vars(contains("num"), aland), sum, na.rm = T) %>%
  distinct() %>%
  mutate(aland_km2 = aland/1000000) %>%
  mutate_at(vars(contains("num")), funs("density" = ./aland_km2)) %>%
  select(matches("density"), MicrotypeID) %>%
  gather(ActivityID, DensityKmSq, num_parks_density:num_jobs_density) %>%
  mutate(ActivityID = str_remove(ActivityID, "num_"), ActivityID = str_remove(ActivityID, "_density")) %>%
  mutate(TripPurposeID = case_when( # create link to trip purposes
    ActivityID %in% c("banks", "creditunion", "parks", "childcare", "snap") ~ 'other',
    ActivityID %in%  c("hosp", "pharm", "urgentcare") ~ 'medical',
    ActivityID == "schools" ~ 'school',
    ActivityID == "jobs" ~ 'work'))

# replace missing values with zero
ops$DensityKmSq[is.na(ops$DensityKmSq) | ops$DensityKmSq == Inf] <- 0
ops$DensityKmSq <- sprintf(ops$DensityKmSq, fmt = '%#.5f')

fwrite(ops, file.path(modeldir, "ActivityDensity.csv"), row.names = F)

##################
# Freight demand by microtype
# FreightDemand.csv
# MicrotypeID Mode VMTPerHour
###########################
freight = row2 %>% 
  right_join(tract_labels, by = "GEOID") %>% # Merge with microtype labels
  select(GEOID, vmt_combi_tract, vmt_single_tract, MicrotypeID) %>%
  group_by(MicrotypeID) %>%
  mutate_at(vars(vmt_combi_tract, vmt_single_tract), sum, na.rm = T) %>%
  select(-GEOID) %>% 
  unique() %>%
  mutate(freight_combo = as.integer(vmt_combi_tract/24), #hourly VMT 
         freight_single = as.integer(vmt_single_tract/24)) %>% 
  gather(Mode, VMTPerHour, freight_combo, freight_single) %>%
  select(MicrotypeID, Mode, VMTPerHour) %>% 
  unique() %>%
  arrange(MicrotypeID) %>%
  mutate(VMTPerHour = as.integer(VMTPerHour))

fwrite(freight, file.path(modeldir, "FreightDemand.csv"), row.names = F)

########################
# MicrotypeAssignment.csv @ CCST level
# variables: OriginMicrotypeID	DestinationMicrotypeID	DistanceBinID	(ThroughMicrotypeID)	Portion
# Apportions total number of trips of each distance to each OD microtype pair
#########################
# create data tables for smaller file sizes
# each row is one segment of a trip from Origin tract (o_geoid) to Destination tract (d_geoid), but also includes
# ccst-level characteristics (fid)
# nhts = as.data.table(nhts) %>% # each row is one trip ID
#   filter(!is.na(o_geoid)) %>% 
#   filter(str_detect(o_geoid, "^06")) %>% #drop trips with origin or destination out of state
#   filter(str_detect(d_geoid, "^06")) %>% 
#   arrange(trip_indx) 

# all = lengths %>% # obs = 3,560,026
#   arrange(trip_id, order) %>% 
#   filter(!is.na(order)) %>% # remove rows that do not have origin and destination info or lengths -> #3,399,497 obs
#   rename(FID = thru_fid, 
#          GEOID = thru_geoid) %>% 
#   right_join(nhts, by = c("o_geoid", "d_geoid")) %>% # but there should NOT be multiple routes for the same NHTS trip... so this looks OKAY
#   unique() %>%
#   arrange(trip_indx) %>%
#   filter(!is.na(o_geoid)) %>%
#   group_by(trip_indx) %>%
#   mutate(trip_length_miles = distance/1609.344) %>%
#   ungroup() %>%
#   left_join(tract_labels %>% select(GEOID, MicrotypeID)) %>% #242,617 observations
#   mutate(DistanceBinID = case_when( # Create distance bins to index across
#     trip_length_miles <=1.3  ~ 'short',
#     trip_length_miles >1.3 &  trip_length_miles <=3 ~ 'medium',
#     trip_length_miles > 3 & trip_length_miles <=8 ~ 'long',
#     TRUE ~ 'xlong'))  %>% # Here is where we aggregate from tract (GEOID) to CCTSM (FID) level
#   select(-o_geoid, -d_geoid, -GEOID, -trpmiles, -trip_id, -thru_length_tract) %>% # remove trip distance from NHTS since we now have lengths data
#   arrange(trip_indx, order) %>% 
#   group_by(trip_indx) %>% # keep only rows with highest sequence number for each thru_fid 
#   mutate(flag_fid = case_when( # collapse adjacent tracts within the same FID (but NOT non-adjacent ones)
#     FID == lead(FID) ~ 1, 
#     TRUE ~ 0)) %>%
#   filter(flag_fid == 0) %>%
#   group_by(trip_indx) %>%
#   mutate(order = row_number(),  #re-number sequence so that it is continuous from 1
#          origin = case_when(
#            order == 1 ~ 1,
#            TRUE~ 0),
#          destination  = case_when(
#            order == max(order) ~ 1,
#            TRUE ~ 0)) %>%
#   ungroup() %>%
#   arrange(trip_indx, order) %>%  
#   group_by(trip_indx) %>%
#   mutate(next_geomicrotype = lead(MicrotypeID), # add leading microtype, next one in the sequence
#          OriginMicrotypeID = case_when( # Columns indexing origin and dest microtypes
#            order ==1 ~ `MicrotypeID`),
#          DestMicrotypeID = case_when(
#            destination ==1 ~ `MicrotypeID`)) %>%
#   fill(OriginMicrotypeID, .direction = "down") %>% # Assign same OD for all trips with same trip ID
#   fill(DestMicrotypeID, .direction = "up") %>%
#   ungroup() %>%
#   mutate(thru_length_miles_ccst = thru_length_ccst/1609.344)  %>% # Convert from meters to miles
#   select(-thru_length_ccst, -trip_length_miles, - origin, -destination)  %>% # don't get rid of FID or ORDER
#   mutate(next_geomicrotype = as.character(next_geomicrotype), 
#          next_geomicrotype = case_when(
#            is.na(next_geomicrotype) ~ 'End', # Assign ending microtype
#            TRUE ~ next_geomicrotype  ))

# all = lengths %>% # obs = 3,560,026
#   arrange(trip_id, order) %>% 
#   filter(!is.na(order)) %>% # remove rows that do not have origin and destination info or lengths -> #3,399,497 obs
#   rename(FID = thru_fid, 
#          GEOID = thru_geoid) %>% 
#   right_join(nhts, by = c("o_geoid", "d_geoid")) %>% # but there should NOT be multiple routes for the same NHTS trip... so this looks OKAY
#   unique() %>%
#   arrange(trip_indx) %>%
#   filter(!is.na(o_geoid)) %>%
#   group_by(trip_indx) %>%
#   mutate(trip_length_miles = distance/1609.344) %>%
#   ungroup() %>%
#   left_join(tract_labels %>% select(GEOID, MicrotypeID)) %>% #242,617 observations
#   mutate(DistanceBinID = case_when( # Create distance bins to index across
#     trip_length_miles <=1  ~ 'bin1',
#     trip_length_miles >1 &  trip_length_miles <=2 ~ 'bin2',
#     trip_length_miles >2 &  trip_length_miles <=4 ~ 'bin3',
#     trip_length_miles >4 &  trip_length_miles <=8 ~ 'bin4',
#     trip_length_miles >8 &  trip_length_miles <=15 ~ 'bin5',
#     trip_length_miles >15 &  trip_length_miles <=20 ~ 'bin6',
#     trip_length_miles > 20 & trip_length_miles <=35 ~ 'bin7',
#     TRUE ~ 'bin8'))  %>% # Here is where we aggregate from tract (GEOID) to CCTSM (FID) level
#   select(-o_geoid, -d_geoid, -GEOID, -trpmiles, -trip_id, -thru_length_tract) %>% # remove trip distance from NHTS since we now have lengths data
#   arrange(trip_indx, order) %>% 
#   group_by(trip_indx) %>% # keep only rows with highest sequence number for each thru_fid 
#   mutate(flag_fid = case_when( # collapse adjacent tracts within the same FID (but NOT non-adjacent ones)
#     FID == lead(FID) ~ 1, 
#     TRUE ~ 0)) %>%
#   filter(flag_fid == 0) %>%
#   group_by(trip_indx) %>%
#   mutate(order = row_number(),  #re-number sequence so that it is continuous from 1
#          origin = case_when(
#            order == 1 ~ 1,
#            TRUE~ 0),
#          destination  = case_when(
#            order == max(order) ~ 1,
#            TRUE ~ 0)) %>%
#   ungroup() %>%
#   arrange(trip_indx, order) %>%  
#   group_by(trip_indx) %>%
#   mutate(next_geomicrotype = lead(MicrotypeID), # add leading microtype, next one in the sequence
#          OriginMicrotypeID = case_when( # Columns indexing origin and dest microtypes
#            order ==1 ~ `MicrotypeID`),
#          DestMicrotypeID = case_when(
#            destination ==1 ~ `MicrotypeID`)) %>%
#   fill(OriginMicrotypeID, .direction = "down") %>% # Assign same OD for all trips with same trip ID
#   fill(DestMicrotypeID, .direction = "up") %>%
#   ungroup() %>%
#   mutate(thru_length_miles_ccst = thru_length_ccst/1609.344)  %>% # Convert from meters to miles
#   select(-thru_length_ccst, -trip_length_miles, - origin, -destination)  %>% # don't get rid of FID or ORDER
#   mutate(next_geomicrotype = as.character(next_geomicrotype), 
#          next_geomicrotype = case_when(
#            is.na(next_geomicrotype) ~ 'End', # Assign ending microtype
#            TRUE ~ next_geomicrotype  ))
######################
# MicrotypeAssignment.csv 
# Estimates for each OD-Distance Bin, the proportion of total travel in each through microtype
# Columns = From, To, DistanceBin, ThroughMicrotype, Portion
################

# assign = all %>%
#   rename(From = MicrotypeID) %>%
#   select(-next_geomicrotype, -o_fid, -d_fid, -flag_fid) %>%
#   # now drop the observations in the same trip and FID so we don't double count the thru distance 
#   # now that the same FID appears multiple times in one trip
#   group_by(trip_indx, FID) %>%  
#   filter(order == max(order)) %>%
#   rename(FromMicrotype = OriginMicrotypeID,
#          ToMicrotype = DestMicrotypeID, 
#          ThroughMicrotype = From) %>%
#   group_by(FromMicrotype, ToMicrotype, DistanceBinID) %>%
#   mutate(DEN = sum(thru_length_miles_ccst*wtperfin, na.rm = T)) %>%
#   group_by(FromMicrotype, ToMicrotype, DistanceBinID, ThroughMicrotype) %>%
#   mutate(Portion = sum(wtperfin*thru_length_miles_ccst, na.rm = T)/DEN) %>% 
#   select(FromMicrotype, ToMicrotype, DistanceBinID, ThroughMicrotype, Portion) %>%
#   unique() %>%
#   group_by(FromMicrotype, ToMicrotype, DistanceBinID) %>%
#   mutate(test = sum(Portion)) %>% # Check that portions sum to one
#   select(-test)
# 
# fwrite(assign, file.path(modeldir, "MicrotypeAssignment.csv"), row.names = F)


##################
# AvgTripLengths.csv AT CCST LEVEL
# Variables: MicrotypeID, avg_thru_length
# Export average length spent in each microtype
###################
# same thing here: need to drop the distance spent in the same FID in the same trip
# avg_length = all %>%
#   group_by(trip_indx, FID) %>%  
#   filter(order == max(order)) %>% # drop multiple observations of same FID so we don't double count length
#   ungroup() %>%
#   group_by(MicrotypeID) %>%
#   mutate(avg_thru_length = weighted.mean(thru_length_miles_ccst, w = wtperfin, na.rm = T)) %>%
#   select(MicrotypeID, avg_thru_length) %>%
#   ungroup() %>%
#   unique() %>%
#   filter(!is.na(avg_thru_length))
# 
# fwrite(avg_length, file.path(modeldir, "AvgTripLengths.csv"), row.names = F)

######################
# TRANSITION MATRIX at CCST level
# Generate probability matrix
######################
# Create matrix of only "To" columns
# sample2 = as.data.table(all) %>%
#   select(-distance, - trip_indx, - o_fid, -d_fid, -thru_length_miles_ccst)
# 
# # using data table to create new indicator columns named for the sequence
# # which equal one if the next microtype is equal to i
# geo = c("A", "B", "C", "D", "E", "F")
# 
# for(i in 1:6) for(j in geo) sample2[next_geomicrotype == paste0(j, "_",i), paste0('To',j,'_', i) := 1] # assign value to column
# 
# # Now, count number of occurrences of each sequence, by origin, destination, and trip distance bin
# sample2 = as.data.frame(sample2) %>%
#   rename(From = MicrotypeID) %>%
#   dplyr::mutate(End = case_when( # Add "End" column for number of trip ends
#     next_geomicrotype == "End" & From == DestMicrotypeID ~ 1,
#     TRUE ~ 0))  %>%
#   mutate(across(matches("To"), function(x){sample2$wtperfin*x})) %>% # Add person weights from NHTS, multiplies every indicator by the weight
#   mutate(across(matches("To"), ~replace_na(.x, 0))) %>%
#   group_by(OriginMicrotypeID, DestMicrotypeID, DistanceBinID, From) %>%
#   mutate_at(vars(matches("To"), End), sum) %>% #sum across all thru micro-geo trips for each OD length bin
#   select(-next_geomicrotype) %>%
#   unique() %>%
#   rowwise %>%
#   mutate(DEN = sum(c_across(ToA_1:End)))
# 
# p_matrix = sample2 %>%
#   mutate_at(vars(matches("To")), funs(Prob = ./ DEN)) %>% # Generate portion of total trip segments in each thru microtype
#   select(matches("Prob"), OriginMicrotypeID, DestMicrotypeID, From, DistanceBinID) %>%
#   relocate(OriginMicrotypeID, DestMicrotypeID, DistanceBinID, From) %>% # re-order columns 
#   arrange(OriginMicrotypeID, DestMicrotypeID, DistanceBinID, From) %>%
#   filter(OriginMicrotypeID != 0 & DestMicrotypeID != 0) %>%
#   unique()
# 
# # Remove the "To" prefix and "_Prob" suffix from column names
# names(p_matrix) <- gsub("To","",names(p_matrix))
# names(p_matrix) <- gsub("_Prob","",names(p_matrix))
# names(p_matrix)
# 
# fwrite(p_matrix, file.path(modeldir, "TransitionMatrix.csv"), row.names = F)

