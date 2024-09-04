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

input_dir = 'CleanData_3by5'
output_dir = 'GEMS_inputs_national_3by5'

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
trips = read_csv(os.path.join(input_dir, 'Demand', 'nhts_no_ids_1hrtimebins_with_imputation_3by5.csv'))

pop = read_csv(os.path.join(input_dir, 'Demand', 'nhts_user_classes_inc_veh_sr_3by5.csv'))


#### processed HPMS data ####
network = read_csv(os.path.join(input_dir, 'Network', 'network_microtype_metrics_2.csv'))
network_VMT = read_csv(os.path.join(input_dir, 'Network', 'hpms_vmt_f_system.csv'))

# load bike density input
bike_density = pd.read_csv(os.path.join(input_dir, 'Network', 'bike_availability_2021.csv'))

### externality cost by tract ####
externality = read_csv(os.path.join(input_dir, 'Cost', 'external_costs_mode_tract.csv'))

### highway system cost by tract ####
highway_system_cost = read_csv(os.path.join(input_dir, 'Cost', 'highway_cost_per_tract.csv'))
highway_cost_group = read_csv(os.path.join(input_dir, 'Cost', 'cost_groups_070323.csv'))

### bike system cost by tract ####
bike_system_cost = read_csv(os.path.join(input_dir, 'Cost', 'bike_cost_per_tract.csv'))
### transit system cost by tract ####
transit_system_cost = read_csv(os.path.join(input_dir, 'Cost', 'transit_system_cost.csv'))
print('total transit mileage:')
print(transit_system_cost[['veh_revenue_hours',
'veh_revenue_miles', 'Directional.Route.Miles']].sum())
### user cost by tract ####
driving_user_cost = read_csv(os.path.join(input_dir, 'Cost', 'auto_cost.csv'))
parking_user_cost = read_csv(os.path.join(input_dir, 'Cost', 'parking_tract_2017.csv'))
transit_user_cost = read_csv(os.path.join(input_dir, 'Cost', 'transit_fare_by_tract_2017.csv'))
taxi_user_cost = read_csv(os.path.join(input_dir, 'Cost', 'uber_fare_tract_2017.csv'))

### transit availability input (by dist bin) ###
transit_availability = read_csv(os.path.join(input_dir, 'Demand', 'mode_availability_input_3by5.csv'))

# <codecell>

# Step 1  --  generate 'DistanceBins.csv'  --> 8 bins
# agg level: origin geotype-network microtype, distance bin

dist_bin = [-1, 1, 2, 4, 8, 15, 20, 35, trips['trpmiles'].max()]
dist_bin_label = ['bin1', 'bin2', 'bin3', 'bin4', 'bin5', 'bin6', 'bin7', 'bin8']
dist_bin_mapping = {'dist_under_1':'bin1', 
                  'dist_1-2':'bin2', 
                  'dist_2-4':'bin3', 
                  'dist_4-8':'bin4', 
                  'dist_8-15': 'bin5', 
                  'dist_15-20': 'bin6', 
                  'dist_20-35': 'bin7', 
                  'dist_above_35': 'bin8'}
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

trips.loc[:, 'TripPurposeID'] = trips.loc[:, 'trip_purpose_agg']
# align trip purpose
# trips.loc[trips['TripPurposeID'].isin(["meals", "shopping"]), 'TripPurposeID'] = 'shopping_meals'
# trips.loc[trips['TripPurposeID'].isin(["social"]), 'TripPurposeID'] = 'social'
# trips.loc[trips['TripPurposeID'].isin(["other", "transp_someone"]), 'TripPurposeID'] = 'other'

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


crash_cost_by_typology = crash_cost_by_typology[['geotype', 'network_microtype', 'Mode', 'cost_crash']]
crash_cost_by_typology.to_csv(os.path.join(output_dir, "SafetyExternalities.csv"), index = False)
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
# step 12 -- generate highway network (lane miles only) --> MFD parameters supplied separately
network_with_typology = pd.merge(network, label_2020,
                       left_on = 'tract', right_on = 'GEOID', how = 'left')
network_with_typology = network_with_typology[['GEOID', 'geotype',  'network_microtype','laneMiles']]
sub_network = network_with_typology.groupby(['geotype',  'network_microtype'])[['laneMiles']].sum()
sub_network = sub_network.reset_index()

sub_network.to_csv(os.path.join(output_dir, "SubNetworks.csv"), index = False)
# <codecell>

# step 13 -- generate transit cost
# agg level geotype + network microtypes
# starting transit cost...
# pg.linear_regression(data[['X', 'Z']], data['Y'])

transit_system_cost.loc[:, 'mode_group'] = 'other'
# reference: https://www.transit.dot.gov/ntd/national-transit-database-ntd-glossary#D
# quote: Rail modes (heavy rail (HR), light rail (LR), commuter rail (CR), inclined plane (IP), cable car (CC) and Monorail/Automated guideway (MG))
transit_system_cost.loc[transit_system_cost['Mode'].isin(["MB", "CB", "RB"]), 'mode_group'] = 'bus'
transit_system_cost.loc[transit_system_cost['Mode'].isin(["LR", "CR", "HR", "YR", "SR", "IP", "CC", "MG"]), 'mode_group'] = 'rail'

transit_system_cost_filtered = transit_system_cost.loc[transit_system_cost['mode_group'] != 'other']

transit_system_cost_filtered.loc[:, 'total_cost'] = transit_system_cost_filtered.loc[:, 'total_cap_exp'] + \
    transit_system_cost_filtered.loc[:, 'total_op_exp']

transit_system_cost_filtered = transit_system_cost_filtered.loc[transit_system_cost_filtered['total_cost'] > 0]
# WARNING: THIS IS PER DAY COST
transit_system_cost_filtered.loc[:, 'total_cost'] /= 365

label_2020_geotype = label_2020[['spatial_id', 'geotype']]
label_2020_geotype = label_2020_geotype.drop_duplicates(keep = 'first')
transit_system_cost_with_typology = pd.merge(transit_system_cost_filtered,
                                             label_2020_geotype, 
                                             on = 'spatial_id', how = 'inner')
# WARNING: THIS IS PER DAY OPERATION
transit_system_cost_with_typology.loc[:, 'veh_revenue_hours'] /= 365
transit_system_cost_with_typology.loc[:, 'veh_revenue_miles'] /= 365


# drop NA before regression
transit_system_cost_with_typology.dropna(subset = ['mode_group','veh_revenue_hours', 'total_cost'], inplace = True)
geotype_group_mapping = {'CBSA_1': 'CBSA_1', 'CBSA_2': 'CBSA_2',
                         'NONCBSA_1': 'NONCBSA', 'NONCBSA_2': 'NONCBSA'}

transit_system_cost_with_typology.loc[:, 'geotype_group'] = \
    transit_system_cost_with_typology.loc[:, 'geotype'].map(geotype_group_mapping)
# perform regression
unique_gt = transit_system_cost_with_typology.geotype_group.unique()
sample_transit_cost = transit_system_cost_with_typology.groupby(['geotype', 'mode_group']).size()
print(sample_transit_cost)
lr_mode_results = None
for gt in unique_gt:
    print('processing geotype group = ' + gt)
    # bus cost specific for geotype group and mode
    cost_gt = \
        transit_system_cost_with_typology.loc[transit_system_cost_with_typology['geotype_group'] == gt]
    cost_gt_bus = cost_gt.loc[cost_gt['mode_group'] == 'bus']
    sample_bus = len(cost_gt_bus)
    # rail cost for all geotypes
    cost_gt_rail = transit_system_cost_with_typology.loc[transit_system_cost_with_typology['mode_group'] == 'rail']
    sample_rail = len(cost_gt_rail)
    bus_lr = pg.linear_regression(cost_gt_bus[['veh_revenue_hours']], cost_gt_bus['total_cost'])
    rail_lr = pg.linear_regression(cost_gt_rail[['veh_revenue_hours']], cost_gt_rail['total_cost'])
    
    bus_lr.loc[:, 'geotype_group'] = gt
    bus_lr.loc[:, 'mode'] = 'bus'
    bus_lr.loc[:, 'sample'] = sample_bus
    
    rail_lr.loc[:, 'geotype_group'] = gt
    rail_lr.loc[:, 'mode'] = 'rail'
    rail_lr.loc[:, 'sample'] = sample_rail
    
    lr_mode_results = pd.concat([lr_mode_results, bus_lr, rail_lr])
    # break
lr_mode_results.to_csv(os.path.join(output_dir, "calibration/transit_cost_by_hour_geotype_lm.csv"), index = False)

# <codecell>

# step 14 -- mode specific inputs prep

# mode speed generation
trips.loc[trips['mode'] == 'hv', 'mode'] = 'auto'
trips.loc[trips['mode'] == 'taxi', 'mode'] = 'ridehail'
trips.loc[:, 'weighted_time'] = \
    trips.loc[:, 'wtperfin'] * trips.loc[:, 'inv_time'] * 60 # in sec
trips.loc[:, 'weighted_dist'] = \
    trips.loc[:, 'wtperfin'] * trips.loc[:, 'trpmiles'] * 1609.34 # in meter

mode_speed = \
    trips.groupby(['o_geotype', 'o_network_microtype', 'mode'])[['weighted_dist', 'weighted_time']].sum()
mode_speed.loc[:, 'SpeedInMetersPerSecond'] = \
    mode_speed.loc[:, 'weighted_dist'] / mode_speed.loc[:, 'weighted_time']
# in m/s

mode_speed = mode_speed.reset_index()
mode_speed.drop(columns = ['weighted_dist', 'weighted_time'], inplace = True)
mode_speed.rename(columns = {'o_geotype': 'geotype', 
                             'o_network_microtype': 'network_microtype'}, 
                  inplace = True)

mode_speed.to_csv(os.path.join(output_dir, "calibration/mode_speed.csv"), index = False)


mode_split = \
    trips.groupby(['o_geotype', 'o_network_microtype', 'mode'])[['wtperfin']].sum()
mode_split = mode_split.reset_index()
mode_split = mode_split.loc[mode_split['mode'] != 'ridehail']
mode_split.loc[:, 'fraction'] = \
    mode_split.loc[:, 'wtperfin']/ \
        mode_split.groupby(['o_geotype', 'o_network_microtype'])['wtperfin'].transform('sum')
mode_split.drop(columns = ['wtperfin'], inplace = True)
mode_split.to_csv(os.path.join(output_dir, "calibration/NHTS_mode_split.csv"), index = False) 

trips.loc[:, 'weighted_total_time_hr'] = \
    trips.loc[:, 'wtperfin'] * trips.loc[:, 'inv_time'] / 60 # in hr

trips.loc[:, 'weighted_dist_mi'] = \
    trips.loc[:, 'wtperfin'] * trips.loc[:, 'trpmiles']

trip_travel_time = \
    trips.groupby(['o_geotype', 
                   'o_network_microtype', 
                   'mode'])[['weighted_total_time_hr', 
                             'weighted_dist_mi', 'wtperfin']].sum()
trip_travel_time = trip_travel_time.reset_index()
trip_travel_time.loc[:, 'trip_travel_time (h)'] = \
    trip_travel_time.loc[:, 'weighted_total_time_hr']/trip_travel_time.loc[:, 'wtperfin']
trip_travel_time.loc[:, 'trip_length (mile)'] = \
    trip_travel_time.loc[:, 'weighted_dist_mi']/trip_travel_time.loc[:, 'wtperfin']
trip_travel_time.loc[:, 'avg_speed (mph)'] = \
    trip_travel_time.loc[:, 'weighted_dist_mi']/trip_travel_time.loc[:, 'weighted_total_time_hr']
trip_travel_time.rename(columns = {'wtperfin': 'trip_count'}, inplace = True)
trip_travel_time.to_csv(os.path.join(output_dir, "calibration/NHTS_trip_travel_time.csv"), index = False) 
# <codecell>

# mode user cost generation

# load parking duration
parking_duration = pd.read_csv(os.path.join(input_dir, 'Cost', 'parking_duration.csv'))
parking_duration.rename(columns = {'d_geotype':'geotype',
                                   'd_microtype':'network_microtype'}, inplace = True)

parking_user_cost_with_typology = pd.merge(label_2020, parking_user_cost,
                                           left_on = 'GEOID', right_on = 'tractcode', 
                                            how = 'left')
parking_user_cost_with_typology.loc[:, 'parking'].fillna(0, inplace = True)
parking_cost_agg = \
    parking_user_cost_with_typology.groupby(['geotype', 'network_microtype'])[['parking']].mean()
parking_cost_agg = parking_cost_agg.reset_index()
parking_cost_agg = pd.merge(parking_cost_agg, parking_duration,
                            on = ['geotype', 'network_microtype'], 
                            how = 'left')
parking_cost_agg.loc[:, 'parking'] *= parking_cost_agg.loc[:, 'parking_duration'] 
parking_cost_agg.rename(columns = {'parking': 'PerEndCost'}, inplace = True)

driving_user_cost.drop(columns = ['Unnamed: 0', 'mode', 'total_auto_cost',
        'total_miles'], inplace = True)
driving_user_cost.rename(columns = {'o_geotype': 'geotype', 
                                    'o_microtype': 'network_microtype',  
                                    'CostPerMile': 'PerMileCost'}, inplace = True)

auto_user_cost = pd.merge(driving_user_cost, parking_cost_agg,
                          on = ['geotype', 'network_microtype'],
                          how = 'outer')
auto_user_cost.loc[:, 'mode'] = 'auto'

transit_user_cost.loc[transit_user_cost['mode'].isin(['rail_l', 'rail_c']), 'mode'] = 'rail'
transit_user_cost_with_typology = pd.merge(transit_user_cost, label_2020, 
                                           left_on = 'tractcode', right_on = 'GEOID', 
                                           how = 'left')
transit_cost_agg = \
    transit_user_cost_with_typology.groupby(['geotype', 'network_microtype', 'mode'])[['fare']].mean()
transit_cost_agg = transit_cost_agg.reset_index()
transit_cost_agg.rename(columns = {'fare': 'PerStartCost'}, inplace = True)

ridehail_user_cost_with_typology = pd.merge(taxi_user_cost, label_2020, 
                                           left_on = 'tractcode', right_on = 'GEOID', 
                                           how = 'left')
ridehail_cost_agg = \
    ridehail_user_cost_with_typology.groupby(['geotype', 'network_microtype'])[['Minimum Fare', 'Price per Minute',
           'Price Per Mile']].mean()
ridehail_cost_agg = ridehail_cost_agg.reset_index()
ridehail_cost_agg.rename(columns = {'Minimum Fare': 'PerStartCost',
                                    'Price per Minute': 'PerMinuteCost',
                                    'Price Per Mile': 'PerMileCost'}, 
                         inplace = True)
ridehail_cost_agg.loc[:, 'mode'] = 'ridehail'

user_cost_allmode = \
    pd.concat([auto_user_cost, transit_cost_agg, ridehail_cost_agg])
user_cost_allmode = user_cost_allmode[['geotype', 'network_microtype', 'mode',
                                       'PerStartCost', 'PerMinuteCost', 'PerMileCost', 'PerEndCost']]
user_cost_allmode.to_csv(os.path.join(output_dir, "calibration/mode_cost.csv"), 
                         index = False)

# <codecell>

# step 15 -- writing auto specific inputs
# a. bike input
bike_density.replace([np.inf, -np.inf], np.nan, inplace=True)
bike_density.dropna(inplace = True)
bike_density_with_typology = pd.merge(bike_density, label_2020, 
                                      on = 'GEOID', how = 'left')
bike_density_agg = \
    bike_density_with_typology.groupby(['geotype', 'network_microtype'])[['bike_station_per_ppl_2021']].mean()
bike_density_agg = bike_density_agg.reset_index()
bike_density_agg.rename(columns = {'bike_station_per_ppl_2021': 'BikesPerCapita'}, inplace = True)
bike_density_agg.fillna(0, inplace = True)

bike_attributes = mode_speed.loc[mode_speed['mode'] == 'bike']
bike_attributes.loc[:, 'CapCostPerBike'] = 9000
bike_attributes.loc[:, 'DailyOpCostPerBike'] = 3400/365
bike_attributes.loc[:, 'CapCostPerDock'] = 4500
bike_attributes.loc[:, 'DailyOpCostPerDock'] = 1700/365
bike_attributes.loc[:, 'DocksPerBike'] = 2
bike_attributes.loc[:, 'DocksPerStation'] = 15
bike_attributes.loc[:, 'PerStartCost'] = 0
bike_attributes.loc[:, 'VehicleSize'] = 0.2
bike_attributes.loc[:, 'PerMileCost'] = 0
bike_attributes.loc[:, 'PerEndCost'] = 0
bike_attributes.loc[:, 'DedicatedLanePreference'] = 0.7

bike_attributes = pd.merge(bike_attributes, bike_density_agg,
                           on = ['geotype', 'network_microtype'], how = 'left')
bike_attributes.drop(columns = ['mode'], inplace = True)
bike_attributes.to_csv(os.path.join(output_dir, "mode/bike.csv"), 
                         index = False)

# <codecell>

# b. auto input
auto_attributes = user_cost_allmode.loc[user_cost_allmode['mode'] == 'auto']
auto_attributes.fillna(0, inplace = True)
auto_speed = mode_speed.loc[mode_speed['mode'] == 'auto']
auto_attributes = pd.merge(auto_attributes, auto_speed,
                           on = ['geotype', 'network_microtype', 'mode'],
                           how = 'left')

auto_attributes.drop(columns = ['mode'], inplace = True)
auto_attributes.to_csv(os.path.join(output_dir, "mode/auto.csv"), 
                         index = False)

# <codecell>

# c-d. bus + rail input
transit_attributes = \
    user_cost_allmode.loc[user_cost_allmode['mode'].isin(['bus', 'rail'])]
transit_attributes.fillna(0, inplace = True)
transit_speed = mode_speed.loc[mode_speed['mode'].isin(['bus', 'rail'])]
transit_attributes = pd.merge(transit_speed, transit_attributes, 
                           on = ['geotype', 'network_microtype', 'mode'],
                           how = 'outer')

transit_attributes.loc[:, 'geotype_group'] = \
    transit_attributes.loc[:, 'geotype'].map(geotype_group_mapping)

lr_mode_results_sel = lr_mode_results.loc[lr_mode_results['names'] == 'veh_revenue_hours']
opcost = lr_mode_results_sel[['coef', 'geotype_group', 'mode']]
opcost.rename(columns = {'coef': 'VehicleOperatingCostPerHour'}, inplace = True)
transit_attributes = pd.merge(transit_attributes,
                              opcost, on = ['geotype_group', 'mode'], 
                              how = 'left')

population_fraction = \
    population.groupby(['h_geotype', 'h_network_microtype'])[['Population']].sum()
population_fraction = population_fraction.reset_index()
population_fraction.loc[:, 'pop_fraction'] = \
    population_fraction.loc[:, 'Population']/ \
        population_fraction.groupby('h_geotype')['Population'].transform('sum')

population_fraction.rename(columns = {'h_geotype': 'geotype',
                                      'h_network_microtype': 'network_microtype',
                                      'Population': 'PopulationMicrotypeID'}, inplace = True)

# <codecell>
print('total mileage before assignment:')
print(transit_system_cost_with_typology[['veh_revenue_hours',
'veh_revenue_miles', 'Directional.Route.Miles']].sum())

transit_system_cost_agg = \
    transit_system_cost_with_typology.groupby(['mode_group','geotype'])[['veh_revenue_hours',
    'veh_revenue_miles', 'Directional.Route.Miles']].sum()
transit_system_cost_agg = transit_system_cost_agg.reset_index()
transit_system_cost_agg = pd.merge(population_fraction, transit_system_cost_agg, 
                                   on = 'geotype', how = 'left')

transit_system_cost_agg.loc[:, 'VehicleRevenueMilesPerDay'] = \
    transit_system_cost_agg.loc[:, 'veh_revenue_miles'] * \
        transit_system_cost_agg.loc[:, 'pop_fraction']

transit_system_cost_agg.loc[:, 'VehicleRevenueHoursPerDay'] = \
    transit_system_cost_agg.loc[:, 'veh_revenue_hours'] * \
        transit_system_cost_agg.loc[:, 'pop_fraction']

transit_system_cost_agg.loc[:, 'DirectionalRouteMiles'] = \
    transit_system_cost_agg.loc[:, 'Directional.Route.Miles'] * \
        transit_system_cost_agg.loc[:, 'pop_fraction']
transit_system_cost_agg.rename(columns = {'mode_group': 'mode'}, inplace = True)
transit_system_cost_agg = \
    transit_system_cost_agg[['mode', 'geotype', 'network_microtype',
                             'PopulationMicrotypeID', 'VehicleRevenueMilesPerDay',
                             'VehicleRevenueHoursPerDay', 'DirectionalRouteMiles']]

print('total mileage after assignment:')
print(transit_system_cost_agg[['VehicleRevenueMilesPerDay',
'VehicleRevenueHoursPerDay', 'DirectionalRouteMiles']].sum())
# <codecell>
transit_attributes = pd.merge(transit_system_cost_agg, transit_attributes,
                              on = ['geotype','network_microtype', 'mode'], how = 'outer')

bus_attributes = transit_attributes.loc[transit_attributes['mode'] == 'bus']
bus_attributes.loc[:, 'SeniorFareDiscount'] = 0.5
bus_attributes.loc[:, 'Headway'] = 600
bus_attributes.loc[:, 'VehicleCapacity'] = 40
bus_attributes.loc[:, 'VehicleSize'] = 1
bus_attributes.loc[:, 'InterliningPortion'] = 2
bus_attributes.loc[:, 'CoveragePortion'] = 0.039148384
bus_attributes.loc[:, 'StopSpacing'] = 338.0833333
bus_attributes.loc[:, 'PassengerWait'] = 6.868148516
bus_attributes.loc[:, 'PassengerWaitDedicated'] = 3
bus_attributes.loc[:, 'MinStopTime'] = 15
bus_attributes.loc[:, 'AccessDistanceMultiplier'] = 0.27954271

bus_attributes.drop(columns = ['mode', 'geotype_group'], inplace = True)
bus_attributes.to_csv(os.path.join(output_dir, "mode/bus.csv"), 
                         index = False)


rail_attributes = transit_attributes.loc[transit_attributes['mode'] == 'rail']
rail_attributes.loc[:, 'SeniorFareDiscount'] = 0.5
rail_attributes.loc[:, 'Headway'] = 1800
rail_attributes.loc[:, 'CoveragePortion'] = 0.006287694
rail_attributes.loc[:, 'StopSpacing'] = 1067.833333
rail_attributes.loc[:, 'AcccessDistanceMultiplier'] = 0.377716492999999

rail_attributes.drop(columns = ['mode', 'geotype_group'], inplace = True)
rail_attributes.to_csv(os.path.join(output_dir, "mode/rail.csv"), 
                         index = False)

# <codecell>

# e - walking

walking_attributes = mode_speed.loc[mode_speed['mode'] == 'walk']
walking_attributes.loc[:, 'PerStartCost'] = 0
walking_attributes.loc[:, 'VehicleSize'] = 0.2
walking_attributes.loc[:, 'PerMileCost'] = 0
walking_attributes.loc[:, 'PerEndCost'] = 0

walking_attributes.drop(columns = ['mode'], inplace = True)
walking_attributes.to_csv(os.path.join(output_dir, "mode/walk.csv"), 
                         index = False)

# <codecell>

# b. ridehail input
ridehail_attributes = user_cost_allmode.loc[user_cost_allmode['mode'] == 'ridehail']
ridehail_attributes.fillna(0, inplace = True)
ridehail_speed = mode_speed.loc[mode_speed['mode'] == 'ridehail']
ridehail_attributes = pd.merge(ridehail_attributes, ridehail_speed,
                           on = ['geotype', 'network_microtype', 'mode'],
                           how = 'outer')

ridehail_attributes.drop(columns = ['mode'], inplace = True)
ridehail_attributes.to_csv(os.path.join(output_dir, "mode/ridehail.csv"), 
                         index = False)
