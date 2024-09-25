# Code to calibrate external costs by mode for Transport Geo revised results
  # TRANSIT SYSTEM COSTS
  # EXTERNAL COSTS
  # LANE DEDICATION COSTS
# NATALIE POPOVICH
# BERKELEY NATIONAL LAB
# UPDATED: MAR 8 2021
# Updated Nov 5 2021: estimating bus capital and operating costs separately
# Modified directory path by Hung-Chia Yang on 3/7/2023
####################################
  # Output files from this code:
#DATA 
    # external_costs_mode_microtype.csv
    # TransitService.csv
    # RoadNetworkCosts.csv
#####################################

# set working directory and sub-directories
mywd <- "C:/FHWA/data paper"
setwd(mywd)

#figuredir <- "./Results/Figures"
inputsdir <- "./Output"
datadir <- "./Output"
#tabdir <- './Results/Tables'

#Install/load packages
library(tidyverse)
library(tidyr)
install.packages('stringr')
library(stringr)
install.packages('janitor')
library(janitor)
install.packages('stargazer')
library(stargazer)
install.packages('ggplot2')
library(ggplot2)

#######
# LOAD DATASETS
######### 
# cleaned external cost data 
costs <- read.csv(file.path(datadir,"external_costs_mode_tract.csv")) %>% 
  dplyr::rename(GEOID = tract)

# Transit cost data
a <- read.csv(file.path(datadir,"transit_costs_merged_062723.csv")) %>% 
  mutate(total_cost = total_cap_exp + total_op_exp)  %>%
  select(Mode, spatial_id, total_cost, Mode.VOMS, Directional.Route.Miles, veh_revenue_hours, veh_revenue_miles) %>%
  dplyr::rename(fleet_size = Mode.VOMS) 

# bring in cluster labels and assignment to rural or urban based on FHWA adjusted definitions
# results from the intersection of census tracts with FHWA Urban Adjusted definition
# Cleaned cost group assignment 
inputs <- read.csv(file.path(datadir, "cost_groups_070323.csv")) 

# Transp Geo geotype results to distribute costs across micro-geo pairs
labels <- read.csv(file.path(datadir,"ccst_geoid_key_transp_geo_with_imputation.csv")) %>%
  unite('MicrotypeID', c(geotype,microtype), sep = "_", remove = F)

# ROW lane-miles for weights (LANE DEDICATION)
row <- read.csv(file.path(datadir,"row_no_sample.csv")) %>% 
  mutate_at(vars(contains("lm_tract_fsys")),  funs("percent" = . /lm_all_tract))  

# Cleaned lane dedication and road network costs
road_costs <- read.csv(file.path(datadir, "cleaned_road_costs_070323.csv"))

###########
#  DATA CLEANING
###########
# Transp Geo geotype results to distribute costs across micro-geo pairs
labels <- labels %>%
  merge(inputs, by = "GEOID") %>%
  select(GEOID, MicrotypeID, geotype, microtype, loc_type_label, spatial_id)

# externality cost
out <- costs %>% 
  merge(labels, by = "GEOID") %>%
    filter(loc_type == loc_type_label) %>% # keep only costs associated with correct microtype classification
    select(-loc_type_label)

out <- out %>%
  select(-SafetyScalar, -Mode2, -MicrotypeID, -spatial_id) %>%
  rename(mode_lower_nest = mode)

write.csv(out, file.path(datadir, 'externality_cost_per_mile_by_tract.csv'))
# Assign assign each micro-geo pair to urban or rural
#take average value for each micro-geotype pair, for 54 total values for each cost type and mode
#dist of costs
# t1 <- out %>% 
#   tabyl(MicrotypeID, loc_type) %>% 
#   adorn_percentages("row") %>% 
#   adorn_pct_formatting(digits = 1)
# 
# keep <- out %>% 
#   merge(t1, by = "MicrotypeID") %>%
#   mutate(micro_loc_type = case_when(
#          urban > rural ~ "urban", 
#          TRUE ~ "rural")) %>% # keep only costs associated with microtype location type
#   filter(loc_type == micro_loc_type) %>%
#   select(-rural, -urban, -GEOID, -loc_type) %>% 
#   filter(air_cost_tract != 0 & ghg_cost_tract!= 0 & noise_cost_tract != 0 & crash_cost_tract !=0) %>% # remove zero values
# # aggregate external costs by geo_micro_type and mode 
#   group_by(MicrotypeID, Mode) %>% 
#   mutate(cost_air = mean(air_cost_tract, na.rm = T),
#          cost_ghg = mean(ghg_cost_tract, na.rm = T),
#          cost_noise = mean(noise_cost_tract, na.rm = T),
#          cost_crash = mean(crash_cost_tract, na.rm = T)) %>% 
#   select(MicrotypeID, microtype, geotype, Mode, cost_air, cost_ghg, cost_noise, cost_crash) 
# 
# export <-as.data.frame(distinct(keep))

# assign zero values for modes: RAIL, BIKE, WALK 
# m <- data.frame(matrix(0, ncol = 8, nrow = nrow(export)))
# colnames(m) <- colnames(export)
# m[,1] <- export[,1]
# #m$mode_group <- "bike"
# m$count <- rep(1:3)
# m <- m %>% 
#   mutate(Mode = case_when(
#         count ==2 ~ 'walk',
#         count == 3 ~ 'rail', 
#         TRUE ~ 'bike')) %>%
#   select(-count)
# 
# export = export %>% 
#   rbind(m) %>%
#   filter(geotype != "0")
# 
# write.csv(export, file.path(datadir,"external_costs_mode_microtype.csv"), row.names = F)


##############################
# TRANSIT COSTS CALIBRATION 
##############################
# mode groups
# bus = Bus (MB) + Commuter Bus (CB)
# brt =  Bus Rapid Transit (RB) 
# LR = Light Rail (LR)
# rail = Commuter Rail (CR) + Heavy Rail (HR) + Hybrid Rail (YR) + Streetcar Rail (SR)

#IGNORE: # Ferryboat (FB), Inclined Plane (IP), Monorail/Automated Guideway (MG) Cable Car (CC) Demand Response Taxi (DT)
# Trolleybus (TB), Vanpool (VP), Jitney (JT) ,Aerial Tramway (TR) Publico (PB) Demand Response (DR)

##########
# DATA CLEANING
##########

# transit system costs
system <- labels %>% 
  select(spatial_id, geotype) %>% 
  distinct() %>%
  merge(a, by = "spatial_id") %>% # merge with cluster labels
  filter(Mode %in% c("MB", "CB", "LR", "CR", "HR", "YR", "SR")) %>% # keep only modes used in Task 3 model
  mutate(mode_group = case_when(
    Mode %in% c("MB", "CB") ~ "bus", #Aggregate BTS modes to Task 2 modes
    Mode  == "LR" ~ "light_rail",
    Mode %in% c("CR", "HR", "YR", "SR") ~ "commuter_rail")) %>% 
  distinct %>%
  filter(!is.na(total_cost), 
         fleet_size > 0,
         total_cost > 0) %>%   #drop obs with zero fleet size listed
  mutate(cost_per_car = total_cost/fleet_size) %>%
  distinct()

write.csv(system, file.path(datadir, 'transit_system_cost_per_car_per_year.csv'))
###############################  
# regressions for costs by mode allowed to vary by geotype
# ran these separate by geotype
################################
# bus.df <- filter(system, mode_group == "bus") #& geotype == "F")
# #m.bus <- lm(total_cost/365 ~ fleet_size, 
# #            data = bus.df )
# 
# m.bus <- lm(total_cost/365 ~ geotype + fleet_size + fleet_size*geotype, 
#             data = bus.df )
# stargazer(m.bus, type = 'text', digits = 2)
# 
# m.lr <- lm(total_cost/365 ~ geotype + fleet_size + fleet_size*geotype, 
#            data = subset(system, mode_group == "light_rail"))
# 
# #m.lr <- lm(total_cost/365 ~ fleet_size, 
# #           data = subset(system, mode_group == "light_rail"))
# stargazer(m.lr, type = 'text', digits = 2)
# 
# #m.rail <- lm(total_cost/365 ~ fleet_size, 
# #             data = subset(system, mode_group == "commuter_rail" & geotype %in% c("C", "D", "F")))
# m.rail <- lm(total_cost/365 ~ geotype + fleet_size + fleet_size*geotype, 
#              data = subset(system, mode_group == "commuter_rail"))
# stargazer(m.rail, type = 'text', digits = 2)
# 
# #export table of all mode costs
# stargazer(m.bus, m.rail, m.lr, type = 'latex',
#           out = file.path(datadir,"transit_regression_results.tex"))
# 
# #export for annual report
# stargazer(m.bus, m.rail, m.lr, type = 'text', digits = 2)
# 
# 
# ###################
# # TransitService.csv (input for GEMS)
# # AGGREGRATE DIR. ROUTE MILES, VEH REV HOURS, AND VEH REV MILES PER GEOTYPE
# #########################
# out <- system %>% 
#   group_by(geotype, mode_group) %>%
#   mutate_at(vars(veh_revenue_hours, veh_revenue_miles), sum, na.rm = T) %>%
#   filter(Directional.Route.Miles > 0) %>%
#   group_by(geotype, mode_group) %>%
#   mutate(dir_route_miles = sum(Directional.Route.Miles), na.rm = T) %>%
#   select(geotype, mode_group, veh_revenue_hours, veh_revenue_miles, dir_route_miles ) %>%
#   distinct()
# 
# write.csv(out, file.path(datadir,"TransitService_062723.csv"), row.names = F)


##################
# LANE DEDICATION COSTS
##################

###########
# DATA CLEANING AND MERGING
###########
inputs = inputs %>%
  left_join(select(labels, -loc_type_label), by = "GEOID") %>%
  mutate(cost_group = as.character(cost_group)) 

# Assign each micro-geotype to cost group based on majority of group representation

#Xiaodan's note: this seems very arbitrary
keep <- inputs %>% 
  distinct(GEOID, .keep_all = T) %>% # if tract is in multiple regions, keep only first obs
  group_by(MicrotypeID) %>%
  slice(which.max(table(cost_group))) %>%
  filter(!is.na(geotype))

# Tract assignment to cost groups by microtype 
# estimating homogeneity of assignment to FHWA groups

head(inputs)

# Table: average maximum percentage of tracts in the same cost group
coverage = inputs %>%
  filter(!is.na(geotype)) %>%
  select(cost_group, MicrotypeID) %>%
  group_by(MicrotypeID) %>%
  add_tally() %>%
  group_by(MicrotypeID, cost_group) %>%
  add_tally() %>%
  ungroup() %>%
  mutate(pct = nn/n) %>%
  group_by(MicrotypeID) %>%
  mutate(best = max(pct)) %>%
  select(MicrotypeID, best, n) %>%
  unique() %>%
  ungroup() %>%
  mutate(weight = n/ sum(n), 
         avg_best = weighted.mean(best, weight))
     
############
# PLOTS OF ASSIGNMENT
################

My_Theme = theme(
  axis.text.x = element_text(size = 14),
  axis.text.y = element_text(size = 14),
  axis.title.x = element_text(size = 14),
  axis.title.y = element_text(size = 14),
  legend.text=element_text(size = 16))

# FIGURE 31: Tract assignment to cost groups by microtype 
#png(file.path(figuredir, "microtype_sets.png"), height = 500, width = 700)
#inputs = inputs %>% 
#  mutate(microtype = as.character(microtype))
#ggplot(data = subset(inputs,!is.na(cost_group) & !is.na(geotype)), 
#       aes(microtype, fill = cost_group, na.rm = T)) +
#  geom_bar(position = "fill") +
#  facet_wrap(~geotype) +
#  coord_flip() +
#  labs(x = "Microtype", y = "Percent in set") +
#  scale_fill_brewer(palette = "PuBuGn",
#                    labels = c("Rural Flat", "Rural Rolling", "Rural Mountainous", "Urban 5-50k", "Urban 50-500k", "Urban 500k-1 mil", "Urban > 1 mil")) +
#  theme(legend.position = "bottom", legend.title = element_blank()) #+
#  My_Theme
#dev.off()

# generate weights by proportion of lane miles of each functional system
# calculate average proportion of each functional system in each microtype to generate weights
lm = row %>%
  rename(GEOID = tract) %>%
  merge(inputs, by = "GEOID", all.y = T) %>%
  select(contains("_percent"), GEOID, cost_group, geotype, microtype) %>%
  # select(matches("_mean"), MicrotypeID) %>% 
  # distinct() %>% 
  # #reshape long in order to generate f_sys weights 
  gather(f_sys, fsys_weight, lm_tract_fsys1_percent:lm_tract_fsys7_percent) %>%
  mutate(f_sys = str_remove(f_sys,"lm_tract_fsys")) %>%
  mutate(f_sys = str_remove(f_sys,"_percent"))

lm_with_cost <- lm %>%
  merge(road_costs, by = c("cost_group", "f_sys"), all.x=TRUE) %>%
  mutate(fsys_weight = replace_na(fsys_weight, 0))

lm_with_cost_tract <- lm_with_cost %>% 
  group_by(GEOID) %>%
  mutate(LaneDedicationPerLaneMile = weighted.mean(restripe, w = fsys_weight, na.rm = T)) %>% # same for rural/urban
  mutate(ROWConstructionPerLaneMile = case_when(
    Location == 'RURAL' ~ weighted.mean(add_noobs, w = fsys_weight, na.rm = T), 
    Location == 'URBAN'  ~ weighted.mean(add_obsA, w = fsys_weight, na.rm = T))) %>%
  select(GEOID, LaneDedicationPerLaneMile, ROWConstructionPerLaneMile, cost_group, cost_class, microtype, geotype) %>%
  distinct() %>%
  arrange(GEOID) 

write.csv(lm_with_cost_tract, file.path(datadir, 'roadnetwork_cost_per_mile_whole_lifespan.csv'))
# merge costs to tracts by microtype group
# generate appropriate cost for each microtype based on distribution of lane-miles 
# keep <- keep %>% 
#   merge(lm, by = "MicrotypeID") %>%
#   select(MicrotypeID, matches("fsys"), f_sys, cost_group, microtype) %>%
#   merge(road_costs, by = c("cost_group", "f_sys")) 

# bike costs use only local road costs
# bike <- keep %>% 
#   filter(f_sys == 7) %>% 
#   mutate(Mode = "bike",
#          LaneDedicationPerLaneMile = restripe,
#          ROWConstructionPerLaneMile = case_when(
#            Location == 'RURAL'  ~ add_noobs, 
#            Location == 'URBAN' ~ add_obsA)) %>%
#   select(MicrotypeID, Mode, LaneDedicationPerLaneMile, ROWConstructionPerLaneMile) %>% 
#   distinct()

# Lane Dedication and ROW Expansion Costs
# auto and transit modes use weighted distribution of functional system
# if rural, assume no widening obstacles
# auto <- keep %>% 
#   group_by(MicrotypeID) %>%
#   mutate(LaneDedicationPerLaneMile = weighted.mean(restripe, w = fsys_weight, na.rm = T)) %>% # same for rural/urban
#   mutate(ROWConstructionPerLaneMile = case_when(
#     Location == 'RURAL' ~ weighted.mean(add_noobs, w = fsys_weight, na.rm = T), 
#     Location == 'URBAN'  ~ weighted.mean(add_obsA, w = fsys_weight, na.rm = T))) %>%
#   select(MicrotypeID, LaneDedicationPerLaneMile, ROWConstructionPerLaneMile, cost_group, cost_class, microtype) %>%
#   distinct() %>%
#   arrange(MicrotypeID) %>%
#   # for water tracts, use the most expensive Terrain Restriction for Rural Flat
#   mutate(ROWConstructionPerLaneMile = case_when(
#     is.na(ROWConstructionPerLaneMile) ~ 9748000,
#     TRUE ~ ROWConstructionPerLaneMile), 
#     LaneDedicationPerLaneMile = case_when(
#       is.na(LaneDedicationPerLaneMile) ~1160000,
#       TRUE ~ LaneDedicationPerLaneMile),
#     Mode = "hv") %>% ungroup %>% 
#   select(Mode, MicrotypeID, LaneDedicationPerLaneMile, ROWConstructionPerLaneMile)
# 
# bus <- auto %>% 
#   mutate(Mode = "bus")
# 
# out <-rbind(auto, bus, bike)

write.csv(out, file.path(datadir,"RoadNetworkCosts_070323.csv"), row.names  =F )






