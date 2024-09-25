# Code to remove identifiers from the NHTS data
# Note: this has to be done on the server since it uses the raw data
# NATALIE POPOVICH
# BERKELEY LAB
# JULY 20 2020
# Updated May 24 2021: making time bins 1 hour each
#######################################
# Useful resources:
# Weighting Report: https://nhts.ornl.gov/assets/2017%20NHTS%20Weighting%20Report.pdf
# NHTS 2017 User Guide: https://nhts.ornl.gov/assets/NHTS2017_UsersGuide_04232019_1.pdf

setwd('C:/FHWA/For FHWA folks/trip_generation_and_gems_inputs') # C: drive is on the remote server
library(tidyverse)
library(data.table)

##############
# LOAD DATASETS
################
#NHTS trips with geoids
#load('travelmode.data.national.withgeoid.RData')
#head(TravelMode)

# load transport geography version
#load('C:/Data/NHTS/travelmode.national_transgeo.RData')
#head(TravelMode)

# load new version with new microtypes, Shivam
load('CleanData/NHTS/travelmode.test.data.national_transgeo.RData')
head(TravelMode)

# CCSTs
#ccsts = read.csv('./Results/contiguous_blobs_v3.csv') 

#census data for median incomes
acs <- fread("CleanData/acs_data_tracts.csv") %>% 
  select(med_inc, GEOID) %>%
  mutate(GEOID = str_pad(as.character(GEOID), width = 11, side = "left", pad = "0")) 

# home geoids
hhct <- fread("RawData/NHTS-Tract/hhct.csv") %>%
  mutate(HHCNTYFP = str_pad(HHCNTYFP, width =3, side = "left", pad = "0"),
         HHCT = str_pad(HHCT, width =6, side = "left", pad = "0"),
         GEOID = paste(HHSTFIPS, HHCNTYFP, HHCT, sep = ""),
         GEOID = str_pad(as.character(GEOID), width = 11, side = "left", pad = "0")) 

#bring in tract level population to generate trip rates per person
#NOTE: May 24 2021: we don't use this in the GEMS input generation file so ignore for now
#pop <- fread("./Results/tract_level_cluster_labels.csv") %>% 
# select(GEOID, population, spatial_id)  %>%
#mutate(GEOID = as.character(GEOID, width = 11))

# Read in cluster labels, updated Sep 13 2021 to included imputed clusters
clust = fread("CleanData/ccst_geoid_key_transp_geo_with_imputation.csv")  %>%
  select(GEOID,  geotype, MicrotypeID) %>%   ##g_mcrty THIS DOES NOT EXIST IN the file, shivam note
  mutate(GEOID = str_pad(as.character(GEOID), width = 11, side = "left", pad = "0"),
         microtype = MicrotypeID)

#joining census tract 
tripact <- read.csv("RawData/NHTS-Tract/tripct.csv")
tripact <- unite(tripact, 'o_geoid',ORIG_ST : ORIG_CT, sep='') 
tripact <- unite(tripact, 'd_geoid',DEST_ST : DEST_CT, sep = '') 
tripact <- unite(tripact, trip_indx, HOUSEID:TDTRPNUM, sep = '_',remove = F)
tripact <- select(tripact, -ORIG_COUNTRY, - DEST_COUNTRY)
# TravelMode <- TravelMode %>% 
#   mutate(o_geoid =str_pad(as.character(o_geoid), width = 11, side = "left", pad = "0"), 
#          d_geoid =str_pad(as.character(d_geoid), width = 11, side = "left", pad = "0"))
TravelMode <- right_join(TravelMode, tripact)

######################
# want to generate: avg number of trips per hour at each start time bin
# for each home micro-geotype and each destination geo-microtype pair by trip purpose
#####################################
nhts_od = TravelMode %>% 
  filter(choice == TRUE) %>% 
  select(trip_indx, o_geoid,	d_geoid,	trpmiles,	wtperfin)
write.table(nhts_od, 'CleanData/NHTS/nhts_od_pairs_2017_071023.csv')

# keep only observations where mode was chosen (these are actual trips)
trips = TravelMode %>% 
  filter(choice == TRUE) %>% 
  separate(person_indx, c("HOUSEID", "PERSON_ID")) %>%
  #keep only necessary variables to develop trip rates
  select(-choice, -access_time, - wait_time, -inv_time, -trip_indx)

rm(TravelMode)

# generate hourly time bins
trips = trips %>%
  mutate(start_time_bin = case_when(
         strttime >= 0 & strttime < 100 ~ '0',
         strttime >= 100 & strttime < 200 ~ '1',
         strttime >= 200 & strttime < 300 ~ '2',
         strttime >= 300 & strttime < 400 ~ '3',
         strttime >= 400 & strttime < 500 ~ '4',
         strttime >= 500 & strttime < 600 ~ '5',
         strttime >= 600 & strttime < 700 ~ '6',
         strttime >= 700 & strttime < 800 ~ '7',
         strttime >= 800 & strttime < 900 ~ '8',
         strttime >= 900 & strttime < 1000 ~ '9',
         strttime >= 1000 & strttime < 1100 ~ '10',
         strttime >= 1100 & strttime < 1200 ~ '11',
         strttime >= 1200 & strttime < 1300 ~ '12',
         strttime >= 1300 & strttime < 1400 ~ '13',
         strttime >= 1400 & strttime < 1500 ~ '14',
         strttime >= 1500 & strttime < 1600 ~ '15',
         strttime >= 1600 & strttime < 1700 ~ '16',
         strttime >= 1700 & strttime < 1800 ~ '17',
         strttime >= 1800 & strttime < 1900 ~ '18',
         strttime >= 1900 & strttime < 2000 ~ '19',
         strttime >= 2000 & strttime < 2100 ~ '20',
         strttime >= 2100 & strttime < 2200 ~ '21',
         strttime >= 2200 & strttime < 2300 ~ '22',
         strttime >= 2300 & strttime < 2400 ~ '23'))

#summary stats, number of trip origin tracts that are not represented

# number of trips observed per census tract origin
num_trips <- trips %>%  
  group_by(o_geoid) %>%
  count() %>%  
  mutate(bin = case_when(
    n > 50 ~ 51, 
    TRUE ~ as.numeric(n)))

# total trips by OD-pair
num_trips_od <- trips %>%  
  group_by(o_geoid, d_geoid) %>%
  count() 

num_trips_od_micro <- trips %>%  
  group_by(o_microtype, d_microtype) %>%
  count() 

table(num_trips_od_micro$n)

#52,363 tracts represented by at least one trip
hist(num_trips$bin, breaks = 20)

# number of trips per census tract origin by mode
num_mode <- trips %>% 
  group_by(o_geoid, mode) %>%
  count() %>%
  group_by(o_geoid) %>%
  mutate(unique_modes = n_distinct(mode))

table(num_mode$unique_modes)


# number of trips observed per origin microtype-geotype
num_trips_mg <- trips %>%  
  group_by(o_microtype, o_geotype) %>%
  count() 

# number of trips per origin micro-geotype by mode
num_mode_mg <- trips %>% 
  group_by(o_geotype, o_microtype, mode) %>%
  count() %>%
  group_by(o_geotype, o_microtype) %>%
  mutate(unique_modes = n_distinct(mode)) %>%
  select(o_geotype, o_microtype, unique_modes) %>%
  distinct()

table(num_mode_mg$unique_modes)

# number of trips observed per origin microtype
num_trips_m <- trips %>%  
  group_by(o_microtype) %>%
  count() 

# number of trips per origin micro-geotype by mode
num_mode_m <- trips %>% 
  group_by(o_microtype, mode) %>%
  count() %>%
  group_by(o_microtype) %>%
  mutate(unique_modes = n_distinct(mode))

table(num_mode_m$unique_modes)

num_mode_od_m <- trips %>%  
  group_by(o_microtype, d_microtype, mode) %>%
  count()  %>%
  group_by(o_microtype, d_microtype,) %>%
  mutate(unique_modes = n_distinct(mode)) %>%
  select(o_microtype, d_microtype, unique_modes) %>%
  distinct()

table(num_mode_od_m$unique_modes)

##merge CCSTs with OD pairs to create summary table
# SEP 13 2021: this should be obsolete now since each tract should only be assigned to one CCST
# XIAODAN'S Notes: this part is outdated, I commented out the entire section and don't run this 

# blobs = ccsts %>%
#   mutate(pct_overlap = area_tract/area) %>%
#   # if tracts listed in multiple blobs, assign tract to blog with highest overlap
#   group_by(GEOID) %>%
#   mutate(max_overlap = max(pct_overlap)) %>%
#   filter(max_overlap == pct_overlap | is.na(max_overlap)) %>%
#   mutate(count = row_number()) %>%
#   # for now just keep one observation of each tract, even if it's in multiple blobs 
#   filter(count == 1) %>%
#   select(GEOID, blob_number) %>%
#   mutate( o_geoid = GEOID) 
# 
# trips2 <- trips %>%
#   merge(blobs, by = "o_geoid") %>%
#   mutate(o_blob = blob_number) %>%
#   select(-blob_number)
# 
# blobs$d_geoid <- blobs$GEOID
# 
# trips3 <- trips2 %>%
#   merge(select(blobs, d_geoid, blob_number), by = "d_geoid", all.x = T) %>%
#   mutate(d_blob = blob_number) 
# 
# #number of trips per origin blob
# num_trips_b <- trips3 %>%  
#   group_by(o_blob) %>%
#   count() 
# 
# # number of trips per origin blob by mode
# num_mode_b <- trips3 %>% 
#   group_by(o_blob, mode) %>%
#   count() %>%
#   group_by(o_blob) %>%
#   mutate(unique_modes = n_distinct(mode)) %>%
#   select(o_blob, unique_modes) %>%
#   distinct()
# 
# table(num_mode_b$unique_modes)
# 
# #number of trips by origin dest blob
# num_trips_b_od <- trips3 %>%  
#   group_by(o_blob, d_blob) %>%
#   count() 


###############################
# OPTION 1: TRY USING NHTS PERSON WEIGHTS AS PROVIDED TO AGGREGATE TO REGIONAL AVG TRIP RATES
#############################

#hhct$HHSTFIPS <- str_pad(hhct$HHSTFIPS, width =2, side = "left", pad = "0")
#hhct$HHCNTYFP <- str_pad(hhct$HHCNTYFP, width =3, side = "left", pad = "0")
#hhct$HHCT <- str_pad(hhct$HHCT, width =6, side = "left", pad = "0")
#hhct$GEOID <- paste(hhct$HHSTFIPS, hhct$HHCNTYFP, hhct$HHCT, sep = "")

hhct <- hhct %>% 
  select(GEOID, HOUSEID)# %>%
 # merge(pop , by = "GEOID") # merge population counts with household tracts

#merge trips by home tract
trips = trips %>% 
  merge(hhct, by = "HOUSEID") %>% 
  merge(clust, by = "GEOID") %>%
  mutate(h_geotype = geotype,
         h_microtype = microtype)

#merge with tract-level median income to impute for households missing income
export = trips %>%
  left_join(acs) %>%
  #remove identifiers and export to work with offline
  select(-geotype, -microtype, -GEOID, -o_geoid, -d_geoid, -PERSONID)

# fwrite(export, "./Data/CleanData/nhts_no_ids.csv", row.names = F)
fwrite(export, "CleanData/NHTS/nhts_no_ids_1hrtimebins_with_imputation.csv", row.names = F)

