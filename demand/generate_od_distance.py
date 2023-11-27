# -*- coding: utf-8 -*-
"""
Created on Mon Nov 27 10:13:55 2023

@author: xiaodanxu
"""

import pandas as pd
import geopandas as gpd
import time
import os
from pygris import states
import geopy.distance

path = 'C:/FHWA_R2'
analysis_year = 2020
os.chdir(path)
start_time = time.time()

state_tracts = pd.read_csv('spatial_boundary/CleanData/combined_tracts_' + str(analysis_year) +'.csv')
#state_tracts = pd.DataFrame(gdf.drop(columns='geometry'))
state_tracts.loc[:, 'STATEFP'] = state_tracts.loc[:, 'STATEFP'].astype(str).str.zfill(2)
state_tracts.loc[:, 'GEOID'] = state_tracts.loc[:, 'GEOID'].astype(str).str.zfill(11)
print('finish tract loading')

# <codecell>

list_of_states = state_tracts.STATEFP.unique()
us_states = states(cb = True, year = analysis_year)
state_to_drop = ['AS', 'PR', 'GU', 'VI', 'MP']
us_states = us_states.loc[~us_states['STUSPS'].isin(state_to_drop)]
od_file_path = 'Demand/CleanData/OD'
list_of_od_files = os.listdir(od_file_path)

def get_od_dist(o_lat, o_lon, d_lat, d_lon):
    origin = (o_lat, o_lon)
    dest = (d_lat, d_lon)
    dist = geopy.distance.geodesic(origin, dest)
    return dist

for st in list_of_states:
    st_name = us_states.loc[us_states['STATEFP'] == st, 'STUSPS'].values[0]
    print(st_name)
    od_file =  [x for x in list_of_od_files if st_name in x][0]
    # print(od_file)
    od_data_no_dist = pd.read_csv(os.path.join(od_file_path, od_file))
    od_data_no_dist.loc[:, 'w_tract'] = \
        od_data_no_dist.loc[:, 'w_tract'].astype(str).str.zfill(11)
    od_data_no_dist.loc[:, 'h_tract'] = \
        od_data_no_dist.loc[:, 'h_tract'].astype(str).str.zfill(11)
    state_tracts_work = state_tracts[['GEOID', 'INTPTLAT', 'INTPTLON']]
    state_tracts_work.columns = ['w_tract', 'w_lat', 'w_lon']
    
    state_tracts_home = state_tracts[['GEOID', 'INTPTLAT', 'INTPTLON']]
    state_tracts_home.columns = ['h_tract', 'h_lat', 'h_lon']
    od_data = pd.merge(od_data_no_dist, state_tracts_work,
                               on = 'w_tract', how = 'left')
    od_data = pd.merge(od_data, state_tracts_home,
                               on = 'h_tract', how = 'left')
    od_data.loc[:, 'distance'] = od_data.apply(lambda row : get_od_dist(row['h_lat'], row['h_lon'],
                                  row['h_lat'], row['h_lon']), axis = 1) 
    break