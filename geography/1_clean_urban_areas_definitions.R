# NATALIE POPOVICH
# BERKELEY LAB
# updated  Dec 18, 2023 by Xiaodan Xu

# load packages
library(sf)
library(readxl)
library(dplyr)
library(tidycensus)
library(purrr)
library(tigris)
library(stringr)

# Set working directory
mywd <- "C:/FHWA_R2/spatial_boundary"
setwd(mywd)

datadir <- "RawData"
cleandir <- "CleanData"
analysis_year = 2021
###########
# LOAD DATA
################
# import Urbanized areas and clusters

ua <- get_acs(
  geography = "urban area",
  year = analysis_year,  
  variables = c(population = "B01003_001",
                housing_unit = "B25001_001"),
  output = "wide")
# 2020/2021 has 3592 samples, 2022 only has 2637 samples --> not used

ua <- ua %>% select(GEOID, NAME, populationE, housing_unitE) %>%
  rename(UA_ID = GEOID, UA_NAME = NAME, UA_population = populationE, ua_housing = housing_unitE)

# ua <- st_read(file.path(datadir, "UAS/tl_2019_us_uac10"), 
#               query = "SELECT * FROM tl_2019_us_uac10") 


ua_boundary <- urban_areas(cb = FALSE, year = analysis_year)
plot(st_geometry(ua_boundary))

# load 2020 census tracts
tracts <- read.csv(file.path(cleandir, "combined_tracts_2020.csv"))

#downloaded block /ua crosswalk file from census online: https://www2.census.gov/geo/docs/reference/ua/2020_UA_BLOCKS.txt
# from here: https://www.census.gov/programs-surveys/geography/guidance/geo-areas/urban-rural.html
cbg_ua_crosswalk <- read.csv(file.path('spatial_boundary', datadir, '2020_UA_BLOCKS.txt'), sep = "|")


# Field Name|Field Description
# STATE|Two digit State code
# COUNTY|Three digit County code
# TRACT|Six digit 2020 Census Tract code
# BLOCK|Four digit 2020 Census block code
# GEOID|Fifteen digit geographic identification code made up of the STATE, COUNTY, TRACT and BLOCK codes
# AREALAND|Land area of the 2020 Census block (square meters)
# 2020_HOU|2020 Census housing unit count of the 2020 Census block
# 2020_POP|2020 Census popuation of the 2020 Census block
# 2020_UACE|2020 Census Urban Area code
# 2020_UA_NAME|2020 Census Urban Area name

#############
# DATA CLEANING
##############
cbg_ua_crosswalk <- cbg_ua_crosswalk %>%
  mutate(block_ID = str_pad(GEOID, 15, pad = "0")) %>%
  mutate(GEOID = str_sub(block_ID, start = 1, end=11))

cbg_ua_crosswalk <- cbg_ua_crosswalk %>%
  select(GEOID,X2020_UACE) %>% distinct()

cbg_ua_crosswalk <- cbg_ua_crosswalk %>%
  mutate(X2020_UACE = str_pad(X2020_UACE, 5, pad = "0")) %>%
  rename(UA_ID = X2020_UACE)

tract_ua_df <- merge(cbg_ua_crosswalk, ua, by = c('UA_ID'), all.x = TRUE)

tracts_short <- tracts %>% select(GEOID, NAMELSAD) %>%
  mutate(GEOID = str_pad(GEOID, 11, pad = "0"))

tracts_with_ua_attr <- merge(tracts_short, tract_ua_df, by = 'GEOID', all.x = TRUE)

tracts_with_ua_attr <- tracts_with_ua_attr %>%
  mutate(census_urban_area = ifelse(is.na(UA_ID),0, 1),
         UA_population = ifelse(is.na(UA_population),0, UA_population),
         ua_housing = ifelse(is.na(ua_housing),0, ua_housing) ) %>%
  arrange(GEOID,desc(census_urban_area), desc(UA_population))

# DRIP DUPLICATE --> IF A TRACT BELONGS TO MULTIPLE URBAN AREA, assign it to bigger one
tracts_with_ua_nodup <- tracts_with_ua_attr %>%
  distinct(GEOID, .keep_all = TRUE)


# Create different categories of region types 
tracts_with_ua_nodup = tracts_with_ua_nodup %>% 
  mutate(size_type = case_when(
    UA_population >=1000000 ~ 1,
    UA_population >= 500000 & UA_population < 1000000 ~ 2,
    UA_population >= 200000 & UA_population < 500000 ~ 3, 
    UA_population >= 50000 & UA_population < 200000 ~ 4,
    UA_population >= 5000 & UA_population < 50000 ~ 5, 
    UA_population >= 2500 & UA_population < 5000 ~ 6,
    UA_population < 2500 ~7 )) %>% 
  mutate(fhwa_type = case_when( # FHWA categories for the HERS Final Task 1 report (YEAR)
    UA_population >=1000000 ~ 'major_urbanized',
    UA_population >= 500000 & UA_population < 1000000 ~ 'large_urbanized',
    UA_population >= 50000 & UA_population < 500000 ~ 'small_urbanized',
    UA_population >= 5000 & UA_population < 50000 ~ 'small_urban', 
    UA_population < 5000 ~ 'rural' )) %>% 
  mutate(fhwa_type_num = case_when(
    UA_population >=1000000 ~ 5,
    UA_population >= 500000 & UA_population < 1000000 ~ 4,
    UA_population >= 50000 & UA_population < 500000 ~ 3,
    UA_population >= 5000 & UA_population < 50000 ~ 2, 
    UA_population < 5000  ~ 1 )) %>% 
  mutate(fhwa_urban_area = case_when(
    UA_population >= 50000  ~ 'urbanized_area',
    UA_population >= 5000 & UA_population < 50000 ~ 'small_urban',
    UA_population < 5000 & ua_housing >= 2000 ~ 'other_urban', 
    UA_population < 5000 & ua_housing < 2000 ~ 'rural' ))

table(tracts_with_ua_nodup$fhwa_urban_area, tracts_with_ua_nodup$fhwa_type)


#export urban divisions
write.csv(tracts_with_ua_nodup, file.path(cleandir, paste0("urban_divisions_", analysis_year, ".csv")), row.names = F)
st_write(ua_boundary, file.path(cleandir, paste0("urban_divisions_boundary_", analysis_year, ".geojson")))








