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

os.chdir('C:/FHWA/data paper')

plt.style.use('ggplot')

# load NHTS data with mode choice attributes (no location)
NHTS_trips = pd.read_csv('Output/NHTS_data_with_time_cost.csv')
mode_choice_coeff = pd.read_csv('Output/mode_choice_coefficients.csv')
mode_choice_coeff = mode_choice_coeff.rename(columns ={'Mode': 'mode'})
NHTS_columns = NHTS_trips.columns

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
available_user_class = NHTS_trips['PopulationGroupID'].unique()
NHTS_trips.loc[NHTS_trips['mode'] == 'hv', 'mode'] = 'auto'
NHTS_trips.loc[NHTS_trips['mode'] == 'taxi', 'mode'] = 'ridehail'

NHTS_trips.loc[:, 'PopulationGroupID'] = NHTS_trips.loc[:, 'PopulationGroupID'].map(pop_group_mapping)

# <codecell>

NHTS_trips_with_util = pd.merge(NHTS_trips, mode_choice_coeff,
                                on = ['mode', 'PopulationGroupID'],
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

NHTS_trips_with_util = NHTS_trips_with_util.drop('Unnamed: 0', axis = 1)
NHTS_trips_with_util.to_csv('Output/NHTS_data_with_utility.csv', index = False)

