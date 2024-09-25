# -*- coding: utf-8 -*-
"""
Created on Wed Apr  3 11:08:14 2024

@author: xiaodanxu
"""

import pandas as pd
from pandas import read_csv
import os
from os import listdir
import numpy as np

os.chdir('C:/FHWA_R2')


# load inputs -- the data needs to adopt the same geographic resolution for census tracts
census_year = '2020'

spatial_crosswalk = read_csv('spatial_boundary/CleanData/cleaned_lodes8_crosswalk_with_ID.csv')
network_data = read_csv('Network/CleanData/network_microtype_metrics.csv')
osm_data = read_csv('Network/CleanData/osm_metrics.csv')
crosswalk_2010_2020 = read_csv('spatial_boundary/CleanData/census_tract_crosswalk_2010_2020.csv')

# <codecell>

# keep useful variables only
crosswalk_attr = ['trct', 'cbsa', 'cbsaname']
spatial_crosswalk_sel = spatial_crosswalk[crosswalk_attr]
spatial_crosswalk_sel = spatial_crosswalk_sel.rename(columns = {'trct':'GEOID'})


network_data.loc[:, 'percent_laneMiles_f_system_1_2'] = \
    network_data.loc[:, 'percent_laneMiles_f_system_1'] + network_data.loc[:, 'percent_laneMiles_f_system_1']
network_attr = ['tract','pct_controlp', 
                'percent_laneMiles_f_system_1_2',
                'percent_laneMiles_f_system_3', 
                'percent_laneMiles_f_system_4',
                'percent_laneMiles_f_system_5_7', 
                'lane_mile_density']

network_data_sel = network_data[network_attr]
network_data_sel = network_data_sel.rename(columns = {'tract':'GEOID'})

osm_attr = ['GEOID', 'self_loop_proportion', 'street_density', 'avg_street_length', 'circuity_avg']
osm_data_sel = osm_data[osm_attr]

# tract_cross_attr = ['GEOID_TRACT_20', 'AREALAND_TRACT_20',  'GEOID_TRACT_10', 'AREALAND_TRACT_10']
# crosswalk_2010_2020_sel = crosswalk_2010_2020[tract_cross_attr]

# <codecell>
network_data_sel = network_data_sel.drop_duplicates(subset = ['GEOID'], keep = 'first')
network_output = pd.merge(spatial_crosswalk_sel, network_data_sel, 
                          on = 'GEOID', how = 'left')
network_output = pd.merge(network_output, osm_data_sel, 
                          on = 'GEOID', how = 'left')
# network_output = pd.merge(network_output, crosswalk_2010_2020_sel, 
#                           left_on = 'GEOID', right_on = 'GEOID_TRACT_20', how = 'left')

# <codecell>

# prepare output
# network_output = network_output.drop(columns = ['GEOID_TRACT_20'])
print('The network output has missing values:')
print(network_output.isnull().sum())

network_output.to_csv('Network/CleanData/network_cluster_input.csv', index = False)
