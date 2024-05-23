# Code to clean raw road network and lane dedication costs
# merge with UACE and FHWA population classifications to assign tracts to cost groups

# NATALIE POPOVICH
# BERKELEY NATIONAL LAB
# LAST UPDATED: MARCH 8 2021
####################################

# set working directory and sub-directories
mywd <- "C:/FHWA_R2/Cost"
setwd(mywd)

rawdir <- "./RawData"
datadir <- "./CleanData"

# Install/load packages
# install.packages('openxlsx')
library(openxlsx)
# install.packages('magrittr')
library(magrittr)
# install.packages('dplyr')
library(dplyr)
library(stringr)
library(reshape2)
library(tidyr)
library(ggplot2)
####################
# LOAD DATASETS
####################
# Read in location types from spatial merge of FHWA and census tracts
loc_type <- read.csv(file.path(datadir, "urban_divisions_2021.csv")) %>%
   select(GEOID, census_urban_area, fhwa_type) %>%
   mutate(GEOID = str_pad(as.character(GEOID), 11, pad = "0"))
  # mutate(loc_type_label = "urban") %>%
  # 

# Read in road grade 
grade <- read.csv(file.path(datadir,"network_microtype_metrics_2.csv")) %>%
  rename(GEOID = tract) %>%
  mutate(GEOID = str_pad(as.character(GEOID), 11, pad = "0")) %>%
   select(GEOID, RRS, "laneMiles",                     
          "laneMiles_1.0", "laneMiles_2.0", "laneMiles_3.0",                 
          "laneMiles_4.0", "laneMiles_5.0", "laneMiles_6.0",                 
          "laneMiles_7.0") 
  
grade[is.na(grade)] <- 0 # fill missing

grade <- grade %>%
 merge(loc_type, by = "GEOID", all.x = T) %>%   # merge with location type from FHWA 
  distinct()

################
# FOR URBAN AREAS NEED TO USE UACE CODES TO MATCH FHWA
# read in intersection of tracts and urban areas over 5000 people 
# int <- read.csv(file.path(datadir,"ua_tract_intersect.csv")) %>%
#   select(UACE, GEOID, pop_rgn)

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
    Location == "RURAL" & Category == "Rolling" ~ 'rural_roll',
    Location == "RURAL" & Category == "Mountainous" ~ 'rural_mtn',
    Category == "Large Urbanized" ~ 'urban_large',
    Category == "Major Urbanized" ~ 'urban_major',
    Category == "Small Urbanized" ~ 'urban_small',
    Category == "Small Urban" ~ 'urban_xsmall')) %>%
  mutate(cost_group = case_when( 
    cost_class == 'rural_flat' ~ 1, 
    cost_class == 'rural_roll' ~ 2, 
    cost_class == 'rural_mtn' ~ 3, 
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

write.csv(costs, file.path(datadir, "cleaned_road_costs.csv"), row.names = F)

##################
# Cost categories by terrain, population
#########################

# merge road grade with UACE-tract intersection 
inputs =  grade %>%
  # left_join(int) %>%
  # mutate(loc_type_label = case_when(
  #   is.na(loc_type_label) ~ "rural",
  #   TRUE ~ "urban"
  # )) %>%
  mutate(cost_group= case_when(
    fhwa_type == "rural" & RRS <= 2 ~ 1,  # Rural - Flat (rural from FHWA-adjusted urban definition and grade < 2% ) per HPMS field manual
    fhwa_type == "rural" & RRS > 2 & RRS <=4 ~ 2, # Rural - Rolling (rural and 2 < grade < 5)
    fhwa_type == "rural" & RRS ==5  ~ 3, # Rural - Mountainous (rural and grade > 5)
    fhwa_type == 'small_urban' ~ 4, # Urban - X- Small Urban
    fhwa_type == 'small_urbanized' ~ 5, # Urban - Small Urbanized (50,000 < UA Population < 500,000)
    fhwa_type == 'large_urbanized' ~ 6, # Urban - Large Urbanized (500,000 < UA Population < 1 million)
    fhwa_type == 'major_urbanized' ~7)) %>% # Urban - Major Urbanized (UA Population > 1 million)
  mutate(cost_group = as.character(cost_group))

table(inputs$cost_group)

# Export cleaned cost dataset, which has some tracts assigned to multiple groups 
# if it touches both, selt one
names(inputs)

write.csv(inputs, file.path(datadir, "cost_groups_070323.csv"), row.names = F)


# assign unit cost to tract

network <- inputs %>% select("GEOID", "laneMiles_1.0", "laneMiles_2.0",    
                            "laneMiles_3.0", "laneMiles_4.0", "laneMiles_5.0",
                            "laneMiles_6.0","laneMiles_7.0", "cost_group")

network <- reshape2::melt(network, id=c("GEOID", "cost_group"),
                          variable.name = 'item', value.name = 'lanemiles')
network <- network %>% tidyr::separate(item, c("NAME", "f_sys"))

network <- network %>% select(-NAME)

cost_to_use <- costs %>% select("cost_group", "f_sys", "restripe", "add_obsA",  "add_noobs") 

network_with_cost <- network %>%
  merge(cost_to_use, by = c("cost_group", "f_sys")) %>%
  mutate(restripe = restripe * lanemiles, 
         add_obsA = add_obsA * lanemiles,
         add_noobs = add_noobs * lanemiles) # calculate total cost per construct type

network_with_cost <- network_with_cost %>% 
  group_by(GEOID, cost_group) %>% 
  summarise(lanemiles = sum(lanemiles),
            restripe = sum(restripe),
            add_obsA = sum(add_obsA),
            add_noobs = sum(add_noobs))

network_with_cost <- network_with_cost %>%
  mutate(restripe = restripe / lanemiles, 
         add_obsA = add_obsA / lanemiles,
         add_noobs = add_noobs / lanemiles) # lane mile weighted cost per mile

network_with_cost <- network_with_cost %>%
  mutate(region_type = case_when(
    cost_group == 1 ~ 'Rural - Flat',  # Rural - Flat (rural from FHWA-adjusted urban definition and grade < 2% ) per HPMS field manual
    cost_group == 2 ~ 'Rural - Rolling', # Rural - Rolling (rural and 2 < grade < 5)
    cost_group == 3 ~ 'Rural - Mountainous', # Rural - Mountainous (rural and grade > 5)
    cost_group == 4 ~ 'Small Urban', # Urban - X- Small Urban
    cost_group == 5 ~ 'Small Urbanized', # Urban - Small Urbanized (50,000 < UA Population < 500,000)
    cost_group == 6 ~ 'Large Urbanized', # Urban - Large Urbanized (500,000 < UA Population < 1 million)
    cost_group == 7 ~ 'major urbanized'))

ggplot(network_with_cost, aes(x=region_type, y=restripe)) + 
  geom_boxplot(fill="slateblue", alpha=0.2) + 
  xlab("region type") + ylab('lane dedication') # land dedication

ggplot(network_with_cost, aes(x=region_type, y=add_obsA)) + 
  geom_boxplot(fill="purple", alpha=0.2) + 
  xlab("region type") + ylab('ROW construction (dense area)') # land dedication

ggplot(network_with_cost, aes(x=region_type, y=add_noobs)) + 
  geom_boxplot(fill="orange", alpha=0.2) + 
  xlab("region type") + ylab('ROW construction (other)')# land dedication

write.csv(network_with_cost,  file.path(datadir, "highway_cost_per_tract.csv"), row.names = F)
