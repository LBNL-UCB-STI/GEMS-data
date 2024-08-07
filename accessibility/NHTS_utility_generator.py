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
NHTS_trips = pd.read_csv('output/NHTS_data_with_time_cost.csv')
mode_choice_coeff = pd.read_csv('output/mode_choice_coefficients_anna.csv')
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
NHTS_trips_subset.groupby(['o_geotype','o_network_microtype',
                           'd_geotype', 'd_network_microtype', 
                           'mode', 'mode_available'])[['wtperfin', 'weighted_PMT']].sum()
mode_availability = mode_availability.reset_index()
mode_availability = mode_availability.rename(columns = {'wtperfin': 'weighted_trips'})
mode_availability.to_csv('output/gems/mode_availability_input.csv')

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
# pre-processing and data checking
available_mode = NHTS_trips['mode'].unique()
available_user_class = NHTS_trips['populationgroupid'].unique()
NHTS_trips.loc[NHTS_trips['mode'] == 'hv', 'mode'] = 'auto'
NHTS_trips.loc[NHTS_trips['mode'] == 'taxi', 'mode'] = 'ridehail'

NHTS_trips.loc[:, 'populationgroupid'] = \
    NHTS_trips.loc[:, 'populationgroupid'].map(pop_group_mapping)

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
    NHTS_trips_with_util.loc[:, 'BetaTravelTime'] * NHTS_trips_with_util.loc[:, 'travel_time'] + \
    NHTS_trips_with_util.loc[:, 'BetaMonetaryCost'] * NHTS_trips_with_util.loc[:, 'cost']
    # NHTS_trips_with_util.loc[:, 'BikeShare_Bike'] * NHTS_trips_with_util.loc[:, 'density_pop']

# Xiaodan's note -- this is a temporary drop, those variables will be fixed once we rerun mode choice data prep
# NHTS_trips_with_util = NHTS_trips_with_util.drop(columns = ['strttime', 'start_time_bin'])
NHTS_trips_with_util.to_csv('output/NHTS_data_with_utility.csv', index = False)

# <codecell>
NHTS_trips_with_util = NHTS_trips_with_util.sort_values(by = ['houseid', 'person_id', 'tdtrpnum'])
sample_trip_with_util = NHTS_trips_with_util.head(12000)
sample_trip_with_util.to_csv('output/sample_data_with_utility.csv', index = False)
