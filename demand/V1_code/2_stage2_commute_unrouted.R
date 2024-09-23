##########################################################################
# CODE TO GENERATE EUCLIDEAN DISTANCES OF COMMUTES ACROSS EACH CENSUS TRACT
# NATALIE POPOVICH
# Imports shapefiles of commutes and exports csv files with lengths after intersecting with census tracts

#######################################
# Set working directory
mywd <- "C:/FHWA/For FHWA folks/microtype_input_preparation"
setwd(mywd)

datadir <- "RawData"
cleandir <- "CleanData"

# Load packages
library(sf)
library(data.table)
library(tidycensus)
library(purrr)

################
# LOAD DATA
##############
#c <- st_read("~/Box/FHWA/Data/CleanData/commute_lengths/shapefiles/fips_50_59", 
          #   query = "SELECT * FROM lines_fixed_fips_50_59") 

#trying with ALL the lines at once. see what happens. starting at 3:55pm on Thursday
d <- st_read("./RawData/commute_lengths/shapefiles/fips_30_35", 
             query =  "SELECT * FROM lines_fixed_fips_30_35")

# total length of commute trip
d$trip_length <- st_length(d)    # calculated length (NOTE: THESE SHOULD NOT MATCH THE CALCULATED INTERSECTION LENGTHS)

#NOTES ABOUT SPEED USING SF FUNCTIONS
#I usually use st_join when dealing with point-polygon intersections where geometries don't change. 
#For poly-poly or poly-line etc. it depends on the intended outcome:
#           if the intersected geometries are needed => st_intersection
#           we merely want to know which intersections exist - st_intersects
#The latter is much faster and is used by default in st_join.


##########################
# bring in Census tract spatial boundaries


# set your UC Census API key
# get an API key here: https://api.census.gov/data/key_signup.html
census_api_key("e74b4d8c97989e07245040ac84168a638247af9a", overwrite = TRUE)
options(tigris_use_cache = TRUE)

# get FIPS codes
#us <- unique(fips_codes$state)[1:51]

fips <- c(20:29) # set equal to the same numbers that are being pulled in from the commute lenghtts
# get spatial boundaries for census tracts + population even though we don't need pop
tracts<- reduce(
  map(fips, function(x) {
    get_acs(geography = "tract", variables = "B01003_001", 
            state = x, geometry = TRUE)
  }), 
  rbind
)

#remove extra columns
tracts <- tracts[, -c(2:5)]

#convert sdf into sf object
tracts <- st_as_sf(tracts)

# transform the coords of the smaller into the CRS of the larger object
tracts <- st_transform(tracts , st_crs(d))

# intersect the commute lines with the census data 
int = st_intersection(d, tracts) 

# find out about the length of each line segment
int$thru_length = st_length(int)

#export csv file with lengths of each trip in each tract
fwrite(int, file.path(datadir, "commute_lengths/commutes_20_29.csv"))

###################################
# run the remaining commute lines in batches
######################################
# FIPS CODES 1 - 9
d <- st_read(file.path(datadir, "commute_lengths/shapefiles/fips_1_9"), 
             query =  "SELECT * FROM lines_fixed_fips_1_9")
d$trip_length <- st_length(d)    
fips <- c(1,2,4,5,6,8,9)
tracts<- reduce(
  map(fips, function(x) {
    get_acs(geography = "tract", variables = "B01003_001", 
            state = x, geometry = TRUE)
  }), 
  rbind
)
tracts <- tracts[, -c(2:5)]
tracts <- st_as_sf(tracts)
tracts <- st_transform(tracts , st_crs(d))
int = st_intersection(d, tracts) 
int$thru_length = st_length(int)
write.csv(int, file.path(datadir, "commute_lengths/commutes_1_9.csv"))

# looks like I tried 1_9 but then it didn't work? this file is not in the Box data folder.

rm(int, d, tracts)
# FIPS CODE 30 - 39

d <- st_read(file.path(datadir, "commute_lengths/shapefiles/fips_30_39"), 
             query =  "SELECT * FROM lines_fixed_fips_30_39")
d <- d[,-1]
fips <- c(37)
tracts<- reduce(
  map(fips, function(x) {
    get_acs(geography = "tract", variables = "B01003_001", 
            state = x, geometry = TRUE)
  }), 
  rbind
)
tracts <- tracts[, -c(2:5)]
tracts <- st_as_sf(tracts)
tracts <- st_transform(tracts , st_crs(d))
int = st_intersection(d, tracts) 
int$commute_length = st_length(int)

int <- st_drop_geometry(int)

fwrite(int, file.path(datadir, "commute_lengths/commutes_37.csv"))


# NEW YORK (fips 36)
d <- st_read(file.path(datadir, "commute_lengths/shapefiles/fips_1_9"), 
             query =  "SELECT * FROM lines_fixed_fips_1_9")
d <- d[,-1]
fips <- c(2,4,5,8,9)

tracts<- reduce(
  map(fips, function(x) {
    get_acs(geography = "tract", variables = "B01003_001", 
            state = x, geometry = TRUE)
  }), 
  rbind
)
tracts <- tracts[, -c(2:5)]
tracts <- st_as_sf(tracts)
tracts <- st_transform(tracts , st_crs(d))
int = st_intersection(d, tracts) 
int$thru_length = st_length(int)
int <- st_drop_geometry(int)
fwrite(int, file.path(datadir, "commute_lengths/commutes_2_4_5_8_9.csv"))



# FIPS CODE 9
d <- st_read("~/Box/FHWA/Data/CleanData/commute_lengths/shapefiles/fips_30_39", 
             query =  "SELECT * FROM lines_fixed_fips_30_39")
d <- d[,-1]
fips <- c(30:35)
tracts<- reduce(
  map(fips, function(x) {
    get_acs(geography = "tract", variables = "B01003_001", 
            state = x, geometry = TRUE)
  }), 
  rbind
)
tracts <- tracts[, -c(2:5)]
tracts <- st_as_sf(tracts)
tracts <- st_transform(tracts , st_crs(d))
int = st_intersection(d, tracts) 
int$commute_length = st_length(int)
int <- st_drop_geometry(int)
fwrite(int, file = "./RawData/commute_lengths/commutes_30_35.csv")


rm(d,int,tracts)
d <- st_read("./RawData//commute_lengths/shapefiles/lines_50_56", 
             query =  "SELECT * FROM lines_fixed_fips_50_59")?
d <- d[,-1]
fips <- c(50,51,53:56)
tracts<- reduce(
  map(fips, function(x) {
    get_acs(geography = "tract", variables = "B01003_001", 
            state = x, geometry = TRUE)
  }), 
  rbind
)
tracts <- tracts[, -c(2:5)]
tracts <- st_as_sf(tracts)
tracts <- st_transform(tracts , st_crs(d))
int = st_intersection(d, tracts) 
int$commute_length = st_length(int)
int <- st_drop_geometry(int)
fwrite(int, file = "./RawData/commute_lengths/commutes_50_56.csv")

