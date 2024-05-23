# Code to clean raw road network and lane dedication costs
# merge with UACE and FHWA population classifications to assign tracts to cost groups

# NATALIE POPOVICH
# BERKELEY NATIONAL LAB
# LAST UPDATED: MARCH 8 2021
####################################

# set working directory and sub-directories
mywd <- "C:/FHWA/For FHWA folks/Road_network_and_transit_cost_generation"
setwd(mywd)

rawdir <- "./RawData"
datadir <- "./CleanData"

# Install/load packages
install.packages('openxlsx')
library(openxlsx)
install.packages('magrittr')
library(magrittr)
install.packages('dplyr')
library(dplyr)

####################
# LOAD DATASETS
####################
# Read in location types from spatial merge of FHWA and census tracts
loc_type <- read.csv(file.path(datadir, "tract_fhwa_intersection.csv")) %>%
  select(GEOID, layer) %>%
  mutate(loc_type_label = "urban") %>%
  mutate(GEOID = as.character(GEOID))

# Read in road grade 
grade <- read.csv(file.path(datadir,"microtypes_inputs.csv")) %>%
  rename(GEOID = tract) %>%
  select(GEOID, road_grade) %>% 
 merge(loc_type, by = "GEOID", all.x = T) %>%   # merge with location type from FHWA 
  distinct()

################
# FOR URBAN AREAS NEED TO USE UACE CODES TO MATCH FHWA
# read in intersection of tracts and urban areas over 5000 people 
int <- read.csv(file.path(datadir,"ua_tract_intersect.csv")) %>%
  select(UACE, GEOID, pop_rgn)

#read in RE-STRIPING costs (FOR LANE-DEDICATION COSTS)
l <- openxlsx::read.xlsx(file.path(rawdir,"G_01_AppA_H_TypUrbCapcCostsPerLM_A-8_2018-09-28+.xlsx"), 
                         sheet = "HERS_Highway_Improvement_Costs", colNames = T, skipEmptyRows = T)
l <- set_names(l, nm = l[1, ])
l <- l[-c(1),c(1:4)]

# read in road network capital costs (FOR ROW EXPANSION COSTS)
r <- openxlsx::read.xlsx(file.path(rawdir,"G_01_AppA_H_TypUrbCapcCostsPerLM_A-8_2018-09-28+.xlsx"), 
                         sheet = "HERS_Capacity_Improvement_Costs", colNames = T, skipEmptyRows = T)
r <- set_names(r, nm = r[1, ])
r <- r[-c(1),c(1:5)]

costs <- l %>%
  merge(r, by = c("Location", "Functional Class", "Category"))
rm(l,r)

################
# generate cost sets for raw cost data
###############
costs = costs %>%
  mutate(cost_class= case_when(
    Location == "RURAL" & Category == "Flat" ~ 'rural_flat', 
    Location == "RURAL" & Category == "Mountainous" ~ 'rural_mtn',
    Location == "RURAL" & Category == "Rolling" ~ 'rural_roll',
    Category == "Large Urbanized" ~ 'urban_large',
    Category == "Major Urbanized" ~ 'urban_major',
    Category == "Small Urbanized" ~ 'urban_small',
    Category == "Small Urban" ~ 'urban_xsmall')) %>%
  mutate(cost_group = case_when( 
    cost_class == 'rural_flat' ~ 1, 
    cost_class == 'rural_mtn' ~ 2, 
    cost_class == 'rural_roll' ~ 3, 
    cost_class == 'urban_xsmall' ~ 4, 
    cost_class == 'urban_small' ~ 5,
    cost_class == 'urban_large' ~ 6, 
    cost_class == 'urban_major' ~ 7)) %>%
  mutate(f_sys = case_when( 
    `Functional Class` == 'Interstate' ~ 1, 
    `Functional Class`== 'Other Freeway and Expressway' ~ 2, 
    `Functional Class` == 'Other Principal Arterial' ~ 3, 
    `Functional Class`==  'Minor Arterial' ~ 4, 
    `Functional Class` == 'Major Collector' ~ 5,
    `Functional Class` == 'Minor Collector' ~ 6, 
    `Functional Class` == 'Local' ~ 7)) %>%
  mutate(restripe = as.integer(`Resurface Existing Lane (RESURFACING)`),
         add_obsA = as.integer(`Add Lane, If Obstacle Code A`),
         add_noobs = as.integer(`Add Lane, If No Obstacles (Normal Cost)`))

write.csv(costs, file.path(datadir, "cleaned_road_costs_070323.csv"), row.names = F)

##################
# Cost categories by terrain, population
#########################

# merge road grade with UACE-tract intersection 
inputs =  grade %>%
  left_join(int) %>%
  mutate(loc_type_label = case_when(
    is.na(loc_type_label) ~ "rural",
    TRUE ~ "urban"
  )) %>%
  mutate(cost_group= case_when(
    loc_type_label == "rural" & road_grade < 2 ~ 1,  # Rural - Flat (rural from FHWA-adjusted urban definition and grade < 2% ) per HPMS field manual
    loc_type_label == "rural" & road_grade >= 2 & road_grade <5 ~ 2, # Rural - Rolling (rural and 2 < grade < 5)
    loc_type_label == "rural" & road_grade >=5  ~ 3, # Rural - Mountainous (rural and grade > 5)
    layer == "ua5000" & loc_type_label == "urban" ~ 4, # Urban - X- Small Urban
    pop_rgn > 50000 & pop_rgn < 500000  ~ 5, # Urban - Small Urbanized (50,000 < UA Population < 500,000)
    pop_rgn >= 500000 & pop_rgn < 1000000 ~ 6, # Urban - Large Urbanized (500,000 < UA Population < 1 million)
    pop_rgn >= 1000000 ~7)) %>%
  mutate(cost_group = as.character(cost_group))# Urban - Major Urbanized (UA Population > 1 million)

table(inputs$cost_group)

# Export cleaned cost dataset, which has some tracts assigned to multiple groups 
# if it touches both, selt one
names(inputs)

write.csv(inputs, file.path(datadir, "cost_groups_070323.csv"), row.names = F)
