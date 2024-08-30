# -*- coding: utf-8 -*-
"""
Created on Thu Aug 17 13:27:40 2023

@author: xiaodanxu
"""

# set up python environment
import pandas as pd
import os
from os import listdir
import matplotlib.pyplot as plt
import seaborn as sns
import numpy as np

os.chdir('C:/FHWA_R2/mode_choice_and_demand_generation')

plt.style.use('ggplot')

# load NHTS data with mode choice attributes (no location)
NHTS_trips = pd.read_csv('output/NHTS_data_with_time_cost_3by5.csv')
mode_choice_coeff = pd.read_csv('output/mode_choice_coefficients_implc_v2.csv')
# mode_choice_coeff = mode_choice_coeff.rename(columns ={'Mode': 'mode',
#                                                        'Geotype': 'h_geotype',
#                                                        'IncomeGroup': 'populationgroupid',
#                                                        'TripPurpose': 'trip_purpose_agg'})

NHTS_columns = NHTS_trips.columns

# <codecell>
# calculate mode availability
NHTS_trips_subset = NHTS_trips.loc[NHTS_trips['mode'].isin(['rail', 'bus'])]
NHTS_trips_subset.loc[:, 'weighted_PMT'] = \
    NHTS_trips_subset.loc[:, 'wtperfin'] * NHTS_trips_subset.loc[:, 'trpmiles']
mode_availability = \
pd.pivot_table(NHTS_trips_subset, 
               index = ['houseid', 'person_id', 'tdtrpnum', 'o_geotype','o_network_microtype',
                           'd_geotype', 'd_network_microtype','trip_dist_bin', 'wtperfin'],
               columns = 'mode', values = 'mode_available')
mode_availability = mode_availability.reset_index()
mode_availability.loc[:,'TransitLayer'] = 'bus_rail'
idx_bus_norail = (mode_availability['bus'] == 1) & (mode_availability['rail'] == 0)
idx_nobus_norail = (mode_availability['bus'] == 0) & (mode_availability['rail'] == 0)
idx_nobus_rail = (mode_availability['bus'] == 0) & (mode_availability['rail'] == 1)
mode_availability.loc[idx_bus_norail,'TransitLayer'] = 'bus_norail'
mode_availability.loc[idx_nobus_norail,'TransitLayer'] = 'nobus_norail'
mode_availability.loc[idx_nobus_rail,'TransitLayer'] = 'nobus_rail'

dist_bin_mapping = {'dist_under_1':'bin1', 
                  'dist_1-2':'bin2', 
                  'dist_2-4':'bin3', 
                  'dist_4-8':'bin4', 
                  'dist_8-15': 'bin5', 
                  'dist_15-20': 'bin6', 
                  'dist_20-35': 'bin7', 
                  'dist_above_35': 'bin8'}
mode_availability.loc[:, 'DistanceBinID'] = \
    mode_availability.loc[:, 'trip_dist_bin'].map(dist_bin_mapping)

mode_availability_out = \
    mode_availability.groupby(['o_geotype','o_network_microtype',
                'd_geotype', 'd_network_microtype', 'DistanceBinID', 'TransitLayer'])[['wtperfin']].sum()
mode_availability_out.loc[:, 'Portion'] = \
    mode_availability_out.loc[:, 'wtperfin']/ \
        mode_availability_out.groupby(['o_geotype','o_network_microtype',
                    'd_geotype', 'd_network_microtype', 'DistanceBinID'])['wtperfin'].transform('sum')
mode_availability_out = mode_availability_out.rename(columns = {'wtperfin': 'weighted_trips'})
mode_availability_out.to_csv('output/gems/ModeAvailability.csv')

# <codecell>
pop_group_mapping = {
    'HighIncVehSenior': 'HighIncVeh', 
    'HighIncVeh': 'HighIncVeh', 
    'LowIncVehSenior': 'LowIncVeh', 
    'LowIncVeh': 'LowIncVeh',
    'LowIncNoVeh': 'LowIncNoVeh', 
    'LowIncNoVehSenior': 'LowIncNoVeh', 
    'HighIncNoVeh': 'HighIncNoVeh',
    'HighIncNoVehSenior': 'HighIncNoVeh'
    }

# trip_purp_mapping = {'social': 'social',
#                      'home': 'home',
#                      'work': 'work',
#                      'shopping': 'shopping_meals',
#                      'transp_someone': 'other',
#                      'medical': 'medical',
#                      'meals': 'shopping_meals',
#                      'school': 'school',
#                      'other': 'other'}
# pre-processing and data checking
available_mode = NHTS_trips['mode'].unique()
available_user_class = NHTS_trips['populationgroupid'].unique()
NHTS_trips.loc[NHTS_trips['mode'] == 'hv', 'mode'] = 'auto'
NHTS_trips.loc[NHTS_trips['mode'] == 'taxi', 'mode'] = 'ridehail'

NHTS_trips.loc[:, 'populationgroupid'] = \
    NHTS_trips.loc[:, 'populationgroupid'].map(pop_group_mapping)
    
NHTS_trips.loc[:, 'short_dist_dummy'] = 0
NHTS_trips.loc[:, 'long_dist_dummy'] = 1
#based on mode choice model spec
short_bins = ['dist_under_1',  'dist_1-2']
NHTS_trips.loc[NHTS_trips['trip_dist_bin'].isin(short_bins), 'short_dist_dummy'] = 1
NHTS_trips.loc[NHTS_trips['trip_dist_bin'].isin(short_bins), 'long_dist_dummy'] = 0
# NHTS_trips.loc[:, 'trip_purpose_agg'] = \
#     NHTS_trips.loc[:, 'trip_purpose'].map(trip_purp_mapping)
# <codecell>

NHTS_trips_with_util = pd.merge(NHTS_trips, mode_choice_coeff,
                                on = ['h_geotype', 'mode', 'populationgroupid', 'trip_purpose_agg'],
                                how = 'left')
NHTS_trips_with_util.loc[:, 'travel_time'] = \
    NHTS_trips_with_util.loc[:, 'access_time'] + NHTS_trips_with_util.loc[:, 'wait_time'] \
    + NHTS_trips_with_util.loc[:, 'inv_time']    

var_to_clean = [ 'travel_time', 'cost']
print(len(NHTS_trips_with_util))
NHTS_trips_with_util = NHTS_trips_with_util.dropna(subset = var_to_clean) # no missing values, yay!
print(len(NHTS_trips_with_util))

# apply mode choice utility function
NHTS_trips_with_util.loc[:, 'utility'] = NHTS_trips_with_util.loc[:, 'Intercept'] + \
    NHTS_trips_with_util.loc[:, 'BetaTravelTimeDistBin1Bin2'] * \
        NHTS_trips_with_util.loc[:, 'travel_time'] * NHTS_trips_with_util.loc[:, 'short_dist_dummy'] + \
    NHTS_trips_with_util.loc[:, 'BetaTravelTimeDistBin3Plus'] * \
        NHTS_trips_with_util.loc[:, 'travel_time'] * NHTS_trips_with_util.loc[:, 'long_dist_dummy'] + \
    NHTS_trips_with_util.loc[:, 'BetaMonetaryCostDistBin1Bin2'] * \
        NHTS_trips_with_util.loc[:, 'cost'] * NHTS_trips_with_util.loc[:, 'short_dist_dummy'] + \
    NHTS_trips_with_util.loc[:, 'BetaMonetaryCostDistBin3Plus'] * \
        NHTS_trips_with_util.loc[:, 'cost'] * NHTS_trips_with_util.loc[:, 'long_dist_dummy']        
    # NHTS_trips_with_util.loc[:, 'BikeShare_Bike'] * NHTS_trips_with_util.loc[:, 'density_pop']

# Xiaodan's note -- this is a temporary drop, those variables will be fixed once we rerun mode choice data prep
# NHTS_trips_with_util = NHTS_trips_with_util.drop(columns = ['strttime', 'start_time_bin'])
NHTS_trips_with_util.to_csv('output/NHTS_data_with_utility_3by5.csv', index = False)

# <codecell>
NHTS_trips_with_util = NHTS_trips_with_util.sort_values(by = ['houseid', 'person_id', 'tdtrpnum'])
sample_trip_with_util = NHTS_trips_with_util.head(12000)
sample_trip_with_util.to_csv('output/sample_data_with_utility_3by5.csv', index = False)
