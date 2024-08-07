#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Tue Aug  6 10:13:56 2024

@author: xiaodanxu
"""
import pandas as pd
import os
import numpy as np
from pandas import read_csv
import matplotlib.pyplot as plt
import seaborn as sns
import warnings
import pingouin as pg # for linear regression

warnings.filterwarnings("ignore")

#files exported:

#DONE - DistanceBins.csv  --> 8 bins
#DONE - DistanceDistribution.csv --> 8 bins
#DONE - FreightDemand.csv
#ON SERVER - ModeAvailability.csv
#DONE - OriginDestination.csv 
#DONE - Population.csv 
#AWAITING MODE CHOICE - PopulationGroups.csv
#DONE - TimePeriods.csv 
#ON SERVER - TransitionMatrix.csv
#DONE TripGeneration.csv --> add geotype
#DONE - TripPurposes.csv 

#MODES
#walk.csv
#auto.csv
#rail.csv
#bus.csv
#bike.csv
#ridehail.csv

'''
Part 1 -- load data and perform initial cleaning
'''
data_path = '/Users/xiaodanxu/Library/CloudStorage/GoogleDrive-arielinseu@gmail.com/My Drive/GEMS/accessibility/GEMS_Input'
os.chdir(data_path)

input_dir = 'CleanData'
output_dir = 'GEMS_inputs_national'

#### V2 typology in 2020 census boundary ####
label_2020 = read_csv(os.path.join(input_dir, 'Typology', 'microtype_geotype_output_2020.csv'))

# check unique typology layer combinations 
micro_geotype_combo = label_2020.groupby(['geotype', 'demand_microtype_comb',
       'network_microtype']).size()

microtype_combo = label_2020.groupby(['demand_microtype_comb',
       'network_microtype']).size()

# note: need to fix the variable name from upstream eventually
label_2020 = label_2020.rename(columns = {'demand_microtype_comb': 'demand_microtype'})


#### NHTS data with typology ###
trips = read_csv(os.path.join(input_dir, 'Demand', 'nhts_no_ids_1hrtimebins_with_imputation.csv'))

pop = read_csv(os.path.join(input_dir, 'Demand', 'nhts_user_classes_inc_veh_sr_v2.csv'))


#### processed HPMS data ####
network = read_csv(os.path.join(input_dir, 'Network', 'network_microtype_metrics_2.csv'))
network_VMT = read_csv(os.path.join(input_dir, 'Network', 'hpms_vmt_f_system.csv'))

### externality cost by tract ####
externality = read_csv(os.path.join(input_dir, 'Cost', 'external_costs_mode_tract.csv'))

### highway system cost by tract ####
highway_system_cost = read_csv(os.path.join(input_dir, 'Cost', 'highway_cost_per_tract.csv'))
highway_cost_group = read_csv(os.path.join(input_dir, 'Cost', 'cost_groups_070323.csv'))

### bike system cost by tract ####
bike_system_cost = read_csv(os.path.join(input_dir, 'Cost', 'bike_cost_per_tract.csv'))
### transit system cost by tract ####
transit_system_cost = read_csv(os.path.join(input_dir, 'Cost', 'transit_system_cost.csv'))
# <codecell>

# Step 1  --  generate 'DistanceBins.csv'  --> 8 bins
# agg level: origin geotype-network microtype, distance bin

dist_bin = [-1, 1, 2, 4, 8, 15, 20, 35, trips['trpmiles'].max()]
dist_bin_label = ['bin1', 'bin2', 'bin3', 'bin4', 'bin5', 'bin6', 'bin7', 'bin8']

trips.loc[:, 'DistanceBinID'] = pd.cut(trips['trpmiles'], bins = dist_bin,
                                       labels=dist_bin_label,
                                       right = True)
print(trips.loc[:, 'DistanceBinID'].unique() ) #should have 8 elements

def weighted_average(dataframe, value, weight):
    val = dataframe[value]
    wt = dataframe[weight]
    return (val * wt).sum() / wt.sum()


dist_bin_distribution = \
trips.groupby(['o_geotype', 'o_network_microtype', 'DistanceBinID']).apply(weighted_average, 
                                     'trpmiles', 'wtperfin')

dist_bin_distribution = dist_bin_distribution.reset_index()

# Zach: you can update the variable name here
dist_bin_distribution.columns = ['o_geotype', 'o_network_microtype', 'DistanceBinID', 'MeanDistanceInMiles']

dist_bin_distribution.to_csv(os.path.join(output_dir, "DistanceBins.csv"), index = False)

# <codecell>

# Step 2 -- generate 'Population.csv'
# agg level - home geotype -demand + network microtypes, population group 

population = pop.groupby([ 'h_geotype', 'h_demand_microtype',
                          'h_network_microtype', 'populationgroupid'])[['wtperfin']].sum()
population = population.reset_index()

# Zach: you can update the variable name here
population.columns = ['h_geotype', 'h_demand_microtype', 'h_network_microtype', 
                      'PopulationGroupID', 'Population']

# check unique typology layer combinations  in pop
micro_geotype_pop = population.groupby(['h_geotype', 'h_demand_microtype', 'h_network_microtype']).size()
microtype_pop = population.groupby(['h_demand_microtype', 'h_network_microtype']).size()
population.to_csv(os.path.join(output_dir, "Population.csv"), index = False)

# <codecell>
# Step 3 -- generate TimePeriods.csv
# agg level time bin
time = trips.groupby('start_time_bin').size()
time = time.reset_index()
time.loc[:, 'DurationInHours'] = 1
time = time[['start_time_bin','DurationInHours']]
time.columns = ['TimePeriodID', 'DurationInHours']

time.to_csv(os.path.join(output_dir, "TimePeriods.csv"), index = False)

# <codecell>
# Step 4 -- generate TripPurposes.csv
# agg level: trip purpose

trips.loc[:, 'TripPurposeID'] = trips.loc[:, 'trip_purpose']
# align trip purpose
trips.loc[trips['TripPurposeID'].isin(["meals", "social", "shopping"]), 'TripPurposeID'] = 'leisure'
trips.loc[trips['TripPurposeID'].isin(["other", "transp_someone"]), 'TripPurposeID'] = 'other'

purpose = trips.groupby('TripPurposeID').size()
purpose = purpose.reset_index()
purpose = purpose[['TripPurposeID']]

purpose.to_csv(os.path.join(output_dir, "TripPurposes.csv"), index = False)

# <codecell>

# Step 5 -- generate TripGeneration.csv
# agg level: HomeGeotype, (added) demand microtype, TimePeriodID,PopulationGroupID,TripPurposeID

grouping_var = ['h_geotype', 'h_demand_microtype', 'start_time_bin',
                'TripPurposeID', 'populationgroupid']
pop_grouping_var = ['h_geotype', 'h_demand_microtype', 'populationgroupid']


trips_by_home_demand = trips.groupby(grouping_var)[['wtperfin']].sum()
trips_by_home_demand.columns = ['trip_count']
trips_by_home_demand = trips_by_home_demand.reset_index()



pop_by_home_demand = pop.groupby(pop_grouping_var)[['wtperfin']].sum()
pop_by_home_demand.columns = ['population']
pop_by_home_demand = pop_by_home_demand.reset_index()

ss_by_home_demand = pop.groupby(pop_grouping_var)[['wtperfin']].count()
ss_by_home_demand.columns = ['person sample size']
ss_by_home_demand = ss_by_home_demand.reset_index()
trip_rate_by_home_demand = \
pd.merge(trips_by_home_demand, pop_by_home_demand,
        on = pop_grouping_var, how = 'left')

trip_rate_by_home_demand = \
pd.merge(trip_rate_by_home_demand, ss_by_home_demand,
        on = pop_grouping_var, how = 'left')

trip_rate_by_home_demand = \
pd.merge(trip_rate_by_home_demand, time,
        left_on = 'start_time_bin',right_on = 'TimePeriodID', how = 'left')

trip_rate_by_home_demand.loc[:, 'TripGenerationRatePerHour'] = \
trip_rate_by_home_demand.loc[:, 'trip_count'] / \
trip_rate_by_home_demand.loc[:, 'population'] / \
    trip_rate_by_home_demand.loc[:, 'DurationInHours']

plt.figure(figsize = (10, 5))
g = sns.boxplot(
    data=trip_rate_by_home_demand,
    x="TimePeriodID", y="TripGenerationRatePerHour", 
    hue="h_demand_microtype", showfliers = False, palette = 'Spectral'
)
plt.show()
trip_rate_output = trip_rate_by_home_demand[['h_geotype', 'h_demand_microtype', 
                                             'TripPurposeID', 'populationgroupid', 
                                             'TimePeriodID', 'TripGenerationRatePerHour']] 
# Zach: you can update the variable name here
trip_rate_output = \
    trip_rate_output.rename(columns = {'populationgroupid': 'PopulationGroupID'})
trip_rate_output.to_csv(os.path.join(output_dir, "TripGeneration.csv"), index = False)


# <codecell>

########### PLACEHOLDER FOR MODE CHOICE INPUT ############
# Step 6 -- generate PopulationGroups.csv

# <codecell>

# Step 7 -- generate OriginDestination.csv
# agg level: Home geotype, home demand_microtype, TimePeriodID,	PopulationGroupID, TripPurposeID,	
# Origin geotype and network MicrotypeID, Destination geotype and network MicrotypeID
# assign proportion of trips by OD pairs among each home type, time bin, pop group and trip purpose
trip_agg_var = ['h_network_microtype', 'h_geotype', 'h_demand_microtype',
 'start_time_bin',  'populationgroupid', 'TripPurposeID',
 'o_geotype', 'o_network_microtype', 'd_geotype', 'd_network_microtype']

od_agg_var = ['h_geotype', 'h_demand_microtype', 'h_network_microtype',
              'start_time_bin',  'populationgroupid', 'TripPurposeID']
origin_destination = trips.groupby(trip_agg_var)[['wtperfin']].sum()
origin_destination = origin_destination.reset_index()
origin_destination.loc[:, 'Portion'] = \
    origin_destination.loc[:, 'wtperfin'] / \
        origin_destination.groupby(od_agg_var)['wtperfin'].transform('sum')
        
origin_destination = \
    origin_destination[['h_network_microtype', 'h_geotype', 'h_demand_microtype',
 'start_time_bin',  'populationgroupid', 'TripPurposeID', 
 'o_geotype', 'o_network_microtype', 'd_geotype', 'd_network_microtype','Portion']]

# Zach: you can update the variable name here    
origin_destination = \
    origin_destination.rename(columns = {'start_time_bin': 'TimePeriodID',  
                                         'populationgroupid': 'PopulationGroupID'})
origin_destination.to_csv(os.path.join(output_dir, "OriginDestination.csv"), index = False)

# <codecell>

# Step 8 -- generate DistanceDistribution.csv

# agg level: Origin geotype and network MicrotypeID, 
# Destination geotype and network MicrotypeID, TripPurposeID, DistanceBinID
dist_agg_var = ['o_geotype', 'o_network_microtype', 'd_geotype', 'd_network_microtype',
 'TripPurposeID', 'DistanceBinID']
dist_frac_var = ['o_geotype', 'o_network_microtype', 'd_geotype', 'd_network_microtype',
 'TripPurposeID']
dist_distribution = trips.groupby(dist_agg_var)[['wtperfin']].sum()
dist_distribution = dist_distribution.reset_index()
dist_distribution = dist_distribution.loc[dist_distribution['wtperfin'] > 0]
dist_distribution.loc[:, 'Portion'] = \
    dist_distribution.loc[:, 'wtperfin'] / \
        dist_distribution.groupby(dist_frac_var)['wtperfin'].transform('sum')

dist_distribution = dist_distribution[['o_geotype', 'o_network_microtype', 'd_geotype', 'd_network_microtype',
  'TripPurposeID', 'DistanceBinID', 'Portion']]

dist_distribution.to_csv(os.path.join(output_dir, "DistanceDistribution.csv"), index = False)
        
# <codecell>

# Step 9 -- generate FreightDemand.csv
# agg level geotype -demand + network microtypes
freight_VMT = pd.merge(network_VMT, label_2020,
                       left_on = 'tract', right_on = 'GEOID', how = 'left')

freight_VMT = freight_VMT.groupby(['geotype', 'demand_microtype',
       'network_microtype'])[['vmt_single_unit', 'vmt_combi']].sum()
freight_VMT = freight_VMT.reset_index()
freight_VMT.loc[:, 'freight_combi'] = freight_VMT.loc[:, 'vmt_combi']/24
freight_VMT.loc[:, 'freight_single'] = freight_VMT.loc[:, 'vmt_single_unit']/24
freight_VMT = pd.melt(freight_VMT, id_vars = ['geotype', 'demand_microtype',
       'network_microtype'], value_vars=['freight_combi', 'freight_single'], var_name = 'Mode',
                      value_name='VMTPerHour')
freight_VMT = freight_VMT.reset_index()
freight_VMT.drop(columns = ['index'], inplace = True)
freight_VMT.to_csv(os.path.join(output_dir, "FreightDemand.csv"), index = False)
# calculate hourly VMT

# <codecell>

# step 10 -- generate Externalities.csv
# agg level geotype + network microtypes

network_VMT.loc[:, 'total_VMT'] = \
    network_VMT.loc[:, ['vmt_single_unit', 'vmt_combi',
       'vmt_ldv']].sum(axis = 1)
network_VMT_sel = network_VMT[['tract', 'total_VMT']]
externality_with_vmt = pd.merge(externality, network_VMT, on = 'tract', how = 'left')
externality_with_vmt = pd.merge(externality_with_vmt, label_2020, 
                                left_on = 'tract', right_on = 'GEOID', how = 'left')
air_cost_by_typology = \
externality_with_vmt.groupby(['geotype', 'network_microtype', 'Mode']).apply(weighted_average, 
                                     'air_cost_tract', 'total_VMT')
air_cost_by_typology = air_cost_by_typology.reset_index()
air_cost_by_typology.columns = ['geotype', 'network_microtype', 'Mode', 'cost_air']

ghg_cost_by_typology = \
externality_with_vmt.groupby(['geotype', 'network_microtype', 'Mode']).apply(weighted_average, 
                                     'ghg_cost_tract', 'total_VMT')
ghg_cost_by_typology = ghg_cost_by_typology.reset_index()
ghg_cost_by_typology.columns = ['geotype', 'network_microtype', 'Mode', 'cost_ghg']

noise_cost_by_typology = \
externality_with_vmt.groupby(['geotype', 'network_microtype', 'Mode']).apply(weighted_average, 
                                     'noise_cost_tract', 'total_VMT')
noise_cost_by_typology = noise_cost_by_typology.reset_index()
noise_cost_by_typology.columns = ['geotype', 'network_microtype', 'Mode', 'cost_noise']

crash_cost_by_typology = \
externality_with_vmt.groupby(['geotype', 'network_microtype', 'Mode']).apply(weighted_average, 
                                     'crash_cost_tract', 'total_VMT')
crash_cost_by_typology = crash_cost_by_typology.reset_index()
crash_cost_by_typology.columns = ['geotype', 'network_microtype', 'Mode', 'cost_crash']

#
externality_by_typology = pd.merge(air_cost_by_typology, noise_cost_by_typology,
                                   on = ['geotype', 'network_microtype', 'Mode'], how = 'inner')
externality_by_typology.loc[:, 'PerMileExtCost'] = \
    externality_by_typology.loc[:, 'cost_air'] + externality_by_typology.loc[:, 'cost_noise']

externality_by_typology = externality_by_typology[['geotype', 'network_microtype', 'Mode', 'PerMileExtCost']]
externality_by_typology.to_csv(os.path.join(output_dir, "Externalities.csv"), index = False)

# <codecell>

# step 11 -- generate RoadNetworkCosts.csv
# agg level geotype + network microtypes

highway_system_cost_with_typology = pd.merge(highway_system_cost,label_2020,
                                             on = 'GEOID', how = 'left')
highway_system_cost_with_typology.loc[:, 'Mode'] = 'auto'
bus_system_cost_with_typology = highway_system_cost_with_typology.copy()
bus_system_cost_with_typology.loc[:, 'Mode'] = 'bus'
bike_system_cost_with_typology = pd.merge(bike_system_cost,label_2020,
                                             on = 'GEOID', how = 'left')
bike_system_cost_with_typology.loc[:, 'Mode'] = 'bike'
system_cost_with_typology = pd.concat([highway_system_cost_with_typology,
                                       bus_system_cost_with_typology,
                                       bike_system_cost_with_typology])

system_cost_with_typology = \
    system_cost_with_typology.loc[system_cost_with_typology['lanemiles'] > 0]

daily_scaling = 365 * 20 # 20 year life span, 365 days
system_cost_with_typology.loc[:, 'ROWConstructionPerLaneMile'] = \
    system_cost_with_typology.loc[:, 'add_noobs'] / daily_scaling

# cost group definition as below (copied from R code)
#   mutate(region_type = case_when(
    # cost_group == 1 ~ 'Rural - Flat',  # Rural - Flat (rural from FHWA-adjusted urban definition and grade < 2% ) per HPMS field manual
    # cost_group == 2 ~ 'Rural - Rolling', # Rural - Rolling (rural and 2 < grade < 5)
    # cost_group == 3 ~ 'Rural - Mountainous', # Rural - Mountainous (rural and grade > 5)
    # cost_group == 4 ~ 'Small Urban', # Urban - X- Small Urban
    # cost_group == 5 ~ 'Small Urbanized', # Urban - Small Urbanized (50,000 < UA Population < 500,000)
    # cost_group == 6 ~ 'Large Urbanized', # Urban - Large Urbanized (500,000 < UA Population < 1 million)
    # cost_group == 7 ~ 'major urbanized'))
urban_idx = (system_cost_with_typology['cost_group']>= 4)
system_cost_with_typology.loc[urban_idx, 'ROWConstructionPerLaneMile'] = \
    system_cost_with_typology.loc[urban_idx, 'add_obsA'] / daily_scaling
    
system_cost_with_typology.loc[:, 'LaneDedicationPerLaneMile'] = \
    system_cost_with_typology.loc[:, 'restripe'] / daily_scaling
    
    
paving_cost_by_typology = \
system_cost_with_typology.groupby(['geotype', 'network_microtype', 'Mode']).apply(weighted_average, 
                                     'LaneDedicationPerLaneMile', 'lanemiles')
paving_cost_by_typology = paving_cost_by_typology.reset_index()
paving_cost_by_typology.columns = ['geotype', 'network_microtype', 'Mode', 'LaneDedicationPerLaneMile']

const_cost_by_typology = \
system_cost_with_typology.groupby(['geotype', 'network_microtype', 'Mode']).apply(weighted_average, 
                                     'ROWConstructionPerLaneMile', 'lanemiles')
const_cost_by_typology = const_cost_by_typology.reset_index()
const_cost_by_typology.columns = ['geotype', 'network_microtype', 'Mode', 'ROWConstructionPerLaneMile']


highway_system_by_typology = pd.merge(paving_cost_by_typology, const_cost_by_typology,
                                   on = ['geotype', 'network_microtype', 'Mode'], how = 'inner')

highway_system_by_typology.to_csv(os.path.join(output_dir, "RoadNetworkCosts.csv"), index = False)
# NOTICE: HIGHWAY COST ALREADY IN PER DAY!!!!

# <codecell>

# step 12 -- generate transit cost
# agg level geotype + network microtypes
# starting transit cost...
# pg.linear_regression(data[['X', 'Z']], data['Y'])

transit_system_cost.loc[:, 'mode_group'] = 'rail'
transit_system_cost.loc[transit_system_cost['Mode'].isin(["MB", "CB"]), 'mode_group'] = 'bus'

transit_system_cost.loc[:, 'total_cost'] = transit_system_cost.loc[:, 'total_cap_exp'] + \
    transit_system_cost.loc[:, 'total_op_exp']

transit_system_cost = transit_system_cost.loc[transit_system_cost['total_cost'] > 0]
# WARNING: THIS IS PER DAY COST
transit_system_cost.loc[:, 'total_cost'] /= 365

label_2020_geotype = label_2020[['spatial_id', 'geotype']]
label_2020_geotype = label_2020_geotype.drop_duplicates(keep = 'first')
transit_system_cost_with_typology = pd.merge(transit_system_cost,
                                             label_2020_geotype, 
                                             on = 'spatial_id', how = 'inner')
# WARNING: THIS IS PER DAY OPERATION
transit_system_cost_with_typology.loc[:, 'veh_revenue_hours'] /= 365
transit_system_cost_with_typology.loc[:, 'veh_revenue_miles'] /= 365
unique_gt = transit_system_cost_with_typology.geotype.unique()

# drop NA before regression
transit_system_cost_with_typology.dropna(subset = ['mode_group','veh_revenue_hours', 'total_cost'], inplace = True)

# perform regression
lr_mode_results = None
for gt in unique_gt:
    print('processing geotype = ' + gt)
    cost_gt = \
        transit_system_cost_with_typology.loc[transit_system_cost_with_typology['geotype'] == gt]
    cost_gt_bus = cost_gt.loc[cost_gt['mode_group'] == 'bus']
    cost_gt_rail = cost_gt.loc[cost_gt['mode_group'] == 'rail']
    bus_lr = pg.linear_regression(cost_gt_bus[['veh_revenue_hours']], cost_gt_bus['total_cost'])
    rail_lr = pg.linear_regression(cost_gt_rail[['veh_revenue_hours']], cost_gt_rail['total_cost'])
    
    bus_lr.loc[:, 'geotype'] = gt
    bus_lr.loc[:, 'mode'] = 'bus'
    
    rail_lr.loc[:, 'geotype'] = gt
    rail_lr.loc[:, 'mode'] = 'rail'
    
    lr_mode_results = pd.concat([lr_mode_results, bus_lr, rail_lr])
    # break
lr_mode_results.to_csv(os.path.join(output_dir, "transit_cost_by_hour_geotype_lm.csv"), index = False)