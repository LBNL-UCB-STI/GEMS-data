# code to estimate total amount of PMT in tract and region
# for each CBSA and county (for rural)
# Takes output from R and QGIS that intersected LEHD home-work location distances with
# census tract boundaries to estimate the length of each trip that goes through each tract

#NOTE MAY 13: check that trip numbers arent duplicated for QGIS output and r output.... 

# Updated April 18 2022: trying to clean up and streamline the code a bit... 
# this should ultimately be replaced with the lengths estimated in R one state at a time....
# set working directory
mywd <- "C:/FHWA/For FHWA folks/microtype_input_preparation"
setwd(mywd)

datadir <- "RawData"
cleandir <- "CleanData"

# Load packages
library(sqldf) # use sql commands for csv files 
library(dplyr)
library(data.table)
library(tictoc)
library(mert)
#steps: 
# 2) aggregate total PMT for each spatial ID (CBSA and each county)
# 3) aggregate total PMT for each tract

##########
# LOAD DATA
############
xwalk <- fread(file.path(datadir, "us_xwalk_tract_2017_withID.csv"), header = T) %>%
  select(fips_st, cty, cbsa, GEOID)

# OD pairs for LEHD home and work
vol = fread(file.path(datadir, "od_pairs_tract_2017.csv")) 

# read in commute lengths for each thru tract
# FIRST IMPORT ALL OF THE QGIS FILES SINCE THESE ARE EASY TO MERGE
fp <- list.files(file.path(datadir, "commute_lengths/euclidean/QGIS"), full.names = TRUE,
                 pattern="*.csv", recursive = FALSE)
raw <-lapply(fp,function(f) {fread(f, header = T)})
lengths <- rbindlist(raw, fill = TRUE)
rm(raw)
#some columns have different names, fill in commute length with other columns when applicable
sum(is.na(lengths$commute_length))
lengths$commute_length <- ifelse(is.na(lengths$commute_length), lengths$length_commute, lengths$commute_length)
lengths$commute_length <- ifelse(is.na(lengths$commute_length), lengths$length_thru, lengths$commute_length)
lengths$commute_length <- ifelse(is.na(lengths$commute_length), lengths$length, lengths$commute_length)
sum(is.na(lengths$commute_length))

lengths = lengths %>%
  select(trip_id, GEOID, commute_length)

# merge with number of trips by trip ID
ft <- list.files(file.path(datadir, "commute_lengths/inputs"), full.names = TRUE,
                 pattern="*.csv", recursive = FALSE)
a =  lapply(ft,function(f) {fread(f, header = T)})

trips = rbindlist(a, fill = T) %>%
  select(trip_id, num_trips) %>% 
  unique()
rm(a)

# assign number of trips to each trip ID even though for this metric we don't care about OD pairs, just total PMT
# merge trip lengths with number of trips
lengths = lengths %>%
  merge(trips, by = "trip_id", allow.cartesian = TRUE)
rm(trips)
###############
# bring in commutes that begin and end in same census tract to assign trip lengths to those trips
###############
vol = vol %>%
  filter(h_tract == w_tract) %>%
  select(-V1) %>%
  rename(tract = w_tract, 
         tract_h = h_tract, 
         num_trips = S000)

###################
# import areas for tracts
################
a <- fread(file.path(datadir, "tract_areas_2017.csv")) %>%
  select(-awater2017)

vol = vol %>%
  merge(a, by = "tract") %>%
  # calculate length as one third the square root of the area of census tract
  mutate(commute_length = (1/3)*sqrt(aland2017)) %>%
  select(tract, num_trips, commute_length) %>%
  rename(GEOID = tract)

#generate trip_ids for intra-tract trips
vol$trip_id <- c(50000000:50070777) #make sure it's well past the total number of other trips

# append to other trips
lengths <- rbind(lengths, vol) %>% 
  unique() %>%  # remove duplicates since each tract is counted twice in merge
  mutate(trip_thru_meters = commute_length * num_trips) # generate total PMT for each tract based on number of trips and commute _length (thru_length)

#################
# AGGREGATING PMT BY REGION AND CENSUS TRACT
##############
new <- lengths %>%
  merge(xwalk, by = "GEOID") %>% # aggregate total PMT (meters) by spatial ID 
  group_by(spatial_id) %>% 
  mutate(pmt_region = sum(trip_thru_meters)/1609) %>% 
  group_by(GEOID) %>% # aggregate PMt by census tract
  mutate(pmt_tract = sum(trip_thru_meters)/1609) 

#export tract level summary stats for second-stage analysis
#export <- new[,c(1:3,8:10)] # what is going on here?
export <- export %>% 
  distinct()

fwrite(export, file.path(cleandir,"centricity_QGIS2.csv"), row.names = FALSE, sep = ",")

###########################
# FOR R OUTPUT WHICH IS SLIGHTLY DIFFERENT
################################
fl <- list.files(file.path(datadir, "commute_lengths/euclidean/R"), full.names = TRUE,
                 pattern="*.csv", recursive = FALSE)
raw <-lapply(fl,function(f) {fread(f, header = T)})
lengths <- rbindlist(raw, fill = TRUE)
rm(raw)

#some columns have different names, fill in commute length with other columns when applicable
lengths$commute_length <- ifelse(is.na(lengths$commute_length), lengths$thru_length, lengths$commute_length)

# only keep necessary variables
lengths = lengths %>%
  select(trip_id, GEOID, commute_length)

# merge with number of trips by trip ID
ft <- list.files(file.path(datadir, "commute_lengths/inputs"), full.names = TRUE,
                 pattern="*.csv", recursive = FALSE)
a <-  lapply(ft,function(f) {fread(f, header = T)})
trips <- rbindlist(a, fill = T) %>%
  select(trip_id, num_trips) %>%
  unique()
rm(a)

# assign number of trips to each trip ID even though for this metric we don't care about OD pairs, just total PMT
# merge trip lengths with number of trips
lengths = lengths %>%
  merge(trips, by = "trip_id", allow.cartesian = TRUE)
rm(trips)
###############
# bring in commutes that begin and end in same census tract to assign trip lengths to those trips
###############
vol<- fread("./RawData/od_pairs_tract_2017.csv")
vol <- filter(vol[,2:4], h_tract == w_tract)
colnames(vol) <- c("tract", "tract_h", "num_trips")


vol = vol %>%
  merge(a[,1:2], by = "tract") %>%
  # calculate length as one third the square root of the area of census tract
  mutate(commute_length = (1/3)*sqrt(aland2017)) %>%
  select(tract, num_trips, commute_length) %>%
  rename(GEOID = tract)

rm(a)
#generate trip_ids for intra-tract trips
vol$trip_id <- c(50000000:50070777) #make sure it's well past the total number of other trips

# append to other trips
lengths <- rbind(lengths, vol) %>% 
  unique() %>% # remove duplicates since each tract is counted twice in merge
  # generate total PMT for each tract based on number of trips and commute _length (thru_length)
  mutate(trip_thru_meters = commute_length * num_trips)

#################
# AGGREGATING PMT BY REGION AND CENSUS TRACT
##############
new = lengths %>%
  merge(xwalk, by = "GEOID") %>% 
  group_by(spatial_id) %>% # aggregate total PMT (meters) by spatial ID 
  mutate(pmt_region = sum(trip_thru_meters)/1609) %>%
  ungroup() %>%
  group_by(GEOID) %>% # aggregating PMt by census tract
  mutate(pmt_tract = sum(trip_thru_meters)/1609)

#export tract level summary stats for second-stage analysis
export <- new %>%
  select(GEOID, spatial_id,pmt_region, pmt_tract) %>%
  unique()

fwrite(export, file = file.path(cleandir, "centricity_r_v2.csv"), row.names = FALSE, sep = ",")

# combine the two datasets so there is one comprehensive measure 
a <- fread(file.path(cleandir, "centricity_r_v2.csv"))
b <- fread(file.path(cleandir, "centricity_QGIS2.csv"))

both <- rbind(a,b)

# keep only maximum value 
both <- both[both[, .I[which.max(pmt_tract)], by=GEOID]$V1]
fwrite(both, file.path(cleandir, "centricity_all_v2.csv"), row.names = F)

