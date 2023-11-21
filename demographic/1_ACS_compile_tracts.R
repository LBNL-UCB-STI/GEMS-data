# Generating dataset for geotypes 
# FHWA project 
# Last Updated: 11-20-2023 by Xiaodan Xu, update package, replace variables for vehicle and race, drop housing variables

######################################################################
setwd("C:/FHWA_R2/Demography")
library(stringr)
library(reshape2)
library(dplyr)
library(tidycensus)
library(purrr)
library(tidyr)
library(data.table)
library(totalcensus)
set_path_to_census("~/Desktop/my_census_data")

# for now, only need population and geometry

# set your UC Census API key when you use this for the first time
# get an API key here: https://api.census.gov/data/key_signup.html
# census_api_key("e74b4d8c97989e07245040ac84168a638247af9a", overwrite = TRUE)
# options(tigris_use_cache = TRUE)
readRenviron("~/.Renviron")
analysis_year = 2021
# 

us <- unique(fips_codes$state)[1:51]
#############
# ACS 5 yr 2021 data
###############

# list of variables
vars <- load_variables(analysis_year, "acs5")
write.csv(vars, file.path('CleanData', paste0('ACS_data_dictionary_', analysis_year, '.csv')))
var_tracts <- vars %>% filter(geography == 'tract')
# downloading all of the census data


############################# 
# population level attributes 
############################# 

# population and gender
population_and_gender_age <- reduce(
  purrr::map(us, function(x) {
    get_acs(
    geography = "tract",
    state = x,
    variables = c(population = "B01001_001",
                  pop_male = "B01001_002",
                  pop_m_under5 = "B01001_003",
                  pop_m_5_9 = "B01001_004",
                  pop_m_10_14 = "B01001_005",
                  pop_m_15_17 = "B01001_006",
                  pop_m_18_19 = "B01001_007",
                  pop_m_20 = "B01001_008",
                  pop_m_21 = "B01001_009",
                  pop_m_22_24 = "B01001_010",
                  pop_m_25_29 = "B01001_011",
                  pop_m_30_34 = "B01001_012",
                  pop_m_35_39 = "B01001_013",
                  pop_m_40_44 = "B01001_014",
                  pop_m_45_49 = "B01001_015",
                  pop_m_50_54 = "B01001_016",
                  pop_m_55_59 = "B01001_017",
                  pop_m_60_61 = "B01001_018",
                  pop_m_62_64 = "B01001_019",
                  pop_m_65_66 = "B01001_020",
                  pop_m_67_69 = "B01001_021",
                  pop_m_70_74 = "B01001_022",
                  pop_m_75_79 = "B01001_023",
                  pop_m_80_84 = "B01001_024",
                  pop_m_over85 = "B01001_025",
                  pop_female = "B01001_017",
                  pop_f_under5 = "B01001_027",
                  pop_f_5_9 = "B01001_028",
                  pop_f_10_14 = "B01001_029",
                  pop_f_15_17 = "B01001_030",
                  pop_f_18_19 = "B01001_031",
                  pop_f_20 = "B01001_032",
                  pop_f_21 = "B01001_033",
                  pop_f_22_24 = "B01001_034",
                  pop_f_25_29 = "B01001_035",
                  pop_f_30_34 = "B01001_036",
                  pop_f_35_39 = "B01001_037",
                  pop_f_40_44 = "B01001_038",
                  pop_f_45_49 = "B01001_039",
                  pop_f_50_54 = "B01001_040",
                  pop_f_55_59 = "B01001_041",
                  pop_f_60_61 = "B01001_042",
                  pop_f_62_64 = "B01001_043",
                  pop_f_65_66 = "B01001_044",
                  pop_f_67_69 = "B01001_045",
                  pop_f_70_74 = "B01001_046",
                  pop_f_75_79 = "B01001_047",
                  pop_f_80_84 = "B01001_048",
                  pop_f_over85 = "B01001_049"
),
    year = analysis_year,
    output = "wide")
  }), 
  rbind
) 

race <- reduce(
  purrr::map(us, function(x) {
    get_acs(
      geography = "tract",
      state = x,
      variables = c(white = "B02001_002",
                    african_american = "B02001_003",
                    indian_native = "B02001_004",
                    asian = "B02001_005",
                    pacific_islander = "B02001_006",
                    other = "B02001_007",
                    two_or_more_races = "B02001_008"),
      year = analysis_year,
      output = "wide")
  }), 
  rbind
) 

# EMPLOYMENT STATUS

#search_tablecontents("acs5", years = 2017, keywords = "employment status", view = TRUE)
emp <- reduce(
  purrr::map(us, function(x) {
    get_acs(
      geography = "tract",
      state = x,
      year = analysis_year,  
      variables = c(pop_16_above  = "B23025_001",
                    emp_civil = "B23025_004",
                    emp_army = "B23025_006",
                    emp_unemp = "B23025_005",
                    emp_no_labor_force = "B23025_007"),
      output = "wide")
  }), 
  rbind
) 


edu <- reduce(
  purrr::map(us, function(x) {
    get_acs(
      geography = "tract",
      state = x,
      year = analysis_year,  
      variables = c(pop_age_25_above = "B15003_001",
                    edu_bs = "B15003_022",
                    edu_ms = "B15003_023",
                    edu_prof = "B15003_024",
                    edu_phd = "B15003_025"),
      output = "wide")
  }), 
  rbind
)

# MODE TO WORK
search_tablecontents("acs5", years = 2017, keywords = "B08301", view = TRUE)
commute_mode <- reduce(
  purrr::map(us, function(x) {
    get_acs(
      geography = "tract",
      state = x,
      year = analysis_year,  
      variables = c(workers_age_16_above = "B08301_001",
                    mode_auto = "B08301_002",
                    mode_transit = "B08301_010",
                    mode_taxi = "B08301_016" , 
                    mode_motocycle = "B08301_017" , 
                    mode_bike = "B08301_018" ,
                    mode_walk = "B08301_019", 
                    mode_other = "B08301_020" , 
                    mode_telecommute = "B08301_021"),
      output = "wide")
  }), 
  rbind
)

############################# 
# household level attributes 
############################# 

# hh count and medium hh income
hh_median_income <- reduce(
  purrr::map(us, function(x) {
  get_acs(
  geography = "tract",
  state = x,
  variables = c(
    households = 'B11001_001',
    med_inc_npv = "B19013_001"),
  year = analysis_year,
  output = "wide")
  }), 
rbind
)

# vehicle ownership

vehicles <- reduce(
  purrr::map(us, function(x) {
    get_acs(
      geography = "tract",
      state = x,
      year = analysis_year, 
      variables = c(veh_0 = "B08201_002",
                    veh_1 = "B08201_003",
                    veh_2 = "B08201_004",
                    veh_3 = "B08201_005",
                    veh_4_above = "B08201_006"), 
      output = "wide")
  }), 
  rbind
) 

# POVERTY STATUS

search_tablecontents("acs5", years = 2017, keywords = "poverty status age", view = TRUE)
poverty <- reduce(
  purrr::map(us, function(x) {
    get_acs(
      geography = "tract",
      state = x,
      year = analysis_year, 
      variables = c(pop_report_poverty = "B17001_001",
                    poverty_below = "B17017_002",
                    poverty_above = "B17017_031"), 
      output = "wide")
  }), 
  rbind
) 

##########################
# housing-level attributes (none at this moment)
##########################


  

##########
# merge all datasets 
##############
acs <- Reduce(function(x,y) merge(x = x, y = y, by = c("GEOID", "NAME"), na.omit = FALSE), 
                   list(population_and_gender_age, race, edu, emp,  commute_mode, #population level
                        hh_median_income, poverty, vehicles)) # household level

###########
# export CSV files
###########

fwrite(acs, file = "CleanData/acs_data_tracts_112023.csv", row.names = F)



########### X-Xu note: legacy code/variable kept for future use ##########

# # RACE
# 
# race <- read_acs5year( year = 2017,  states = states_DC,  
#                        table_contents = c("race_total = B03002_001",
#                                           "race_not_latino_total = B03002_002",
#                                           "race_latino_total = B03002_012") ,
#                        summary_level = "tract" )
# race <- race[,c(1,4:12)]
# 

# 
# 
# ###################
# # MOVE IN DATE
# 
# #search_tablecontents("acs5", years = 2017, keywords = "moved into unit", view = TRUE)
# movedin <- read_acs5year( year = 2017, states = states_DC,  
#                           table_contents = c("movedin_total = B25038_001", 
#                                              "moved_post2015o = B25038_003",
#                                              "moved_2010_14o = B25038_004",
#                                              "moved_2000so_ = B25038_005",
#                                              "moved_1990so = B25038_006",
#                                              "moved_1980so = B25038_007",
#                                              "moved_pre1979o = B25038_008",
#                                              "moved_post2015r = B25038_010",
#                                              "moved_2010_14r = B25038_011",
#                                              "moved_2000sr_ = B25038_012",
#                                              "moved_1990sr = B25038_013",
#                                              "moved_1980sr = B25038_014",
#                                              "moved_pre1979r = B25038_015"), 
#                           summary_level = "tract")
# movedin <- movedin[, - c(2,3, 17:22)]
# 
# #################
# # NUMBER OF UNITS IN STRUCTURE
# 
# search_tablecontents("acs5", years = 2017, keywords = "units in structure", view = TRUE)
# units <- read_acs5year( year = 2017,  states = states_DC,   
#                         table_contents = c("structures_total = B25024_001", 
#                                            "structures_units1d = B25024_002", 
#                                            "structures_units1a = B25024_003", 
#                                            "structures_units2 = B25024_004", 
#                                            "structures_units3_4 = B25024_005", 
#                                            "structures_units5_9 = B25024_006", 
#                                            "structures_units10_19 = B25024_007", 
#                                            "structures_units20_49 = B25024_008", 
#                                            "structures_units50plus = B25024_009", 
#                                            "structures_unitsmobile = B25024_010", 
#                                            "structures_unitsother = B25024_011" ), 
#                         summary_level = "tract")
# units <- units[,-c(2,3,15:20)]
# 
# ##########################
# # YEAR STRUCTURE BUILT 
# 
# #search_tablecontents("acs5", years = 2017, keywords = "B25034", view = TRUE)
# built <- read_acs5year( year = 2017,  states = states_DC,   
#                         table_contents = c("built_total = B25034_001", 
#                                            "built_post2014 = B25034_002", 
#                                            "built_2010_13 = B25034_003", 
#                                            "built_2000_09 = B25034_004", 
#                                            "built_1990_99 = B25034_005", 
#                                            "built_1980_89 = B25034_006", 
#                                            "built_1970_79 = B25034_007", 
#                                            "built_1960_69 = B25034_008", 
#                                            "built_1950_59 = B25034_009", 
#                                            "built_1940_49= B25034_010", 
#                                            "built_pre1939 = B25034_011" ), 
#                         summary_level = "tract")
# built <- built[,-c(2,3,15:20)]
# 
# 
# ##########
# # tenure
# search_tablecontents("acs5", years = 2017, keywords = "B25003", view = TRUE)
# tenure <- read_acs5year( year = 2017,  states = states_DC,  
#                          table_contents = c("tenure_total = B25003_001",
#                                             "tenure_own = B25003_002",
#                                             "tenure_rent = B25003_003") ,
#                          summary_level = "tract" )
# tenure <- tenure[,c(1,4:6)]
# 
# ##################
# # VACANCY 
# 
# #search_tablecontents("acs5", years = 2017, keywords = "vacancy", view = TRUE)
# vacancy <- read_acs5year( year = 2017, states = states_DC,    
#                           table_contents = "vacant_units = B25004_001",
#                           summary_level = "tract" )
# vacancy <- vacancy[,c(1,4)]

# # median housing cost
# search_tablecontents("acs5", years = 2017, keywords = "B25064", view = TRUE)
# 
# hcost <- read_acs5year( year = 2017, states = states_DC,    
#                         table_contents = c("median_rent = B25064_001",
#                                            "median_owner_housing_cost = B25088_001"), 
#                         summary_level = "tract" )
# hcost <- hcost[, c(1,4:5)]





