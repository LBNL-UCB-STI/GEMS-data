# -*- coding: utf-8 -*-
"""
Created on Wed Mar  6 11:11:00 2024

@author: xiaodanxu
"""

from pandas import read_csv
import os
from os import listdir
import pandas as pd
import matplotlib.pyplot as plt
import pygris
import geopandas as gpd

os.chdir('C:/FHWA_R2')
data_path = 'RawData'
output_path = 'CleanData'

# NTAD transit data accessed fEB 2024
# LINK: https://data.bts.gov/Bicycles-and-Pedestrians/Locations-of-Docked-Bikeshare-Stations-by-System-a/7m5x-ubud/about_data
bike_share_file = os.path.join('Network', data_path, 'BTS', 'Locations of Docked Bikeshare Stations by System and Year_20240306.geojson')
national_bike_share = gpd.read_file(bike_share_file)


run_type_desc = {1: 'mode choice', 2: 'GEMS update'}
run_type = 2
if run_type == 1:
    ct_file = 'spatial_boundary/CleanData/combined_tracts_2018.geojson'
    analysis_year = 2017
else:
    ct_file = 'spatial_boundary/CleanData/combined_tracts_2020.geojson'
    analysis_year = 2017
 
national_bike_share_sel = national_bike_share.loc[national_bike_share['year'] == str(analysis_year)]
print(len(national_bike_share_sel))

    
us_census_tract = gpd.read_file(ct_file)
# change projection
data_crs = national_bike_share.crs
us_census_tract = us_census_tract.to_crs(data_crs)

# load population data
population_file = os.path.join('Demography', output_path, 'acs_data_tracts_' + str(analysis_year) + '.csv')
acs_population = read_csv(population_file)
# <codecell>


bike_station_to_tract = national_bike_share_sel.sjoin(us_census_tract, how="left")
bike_station_to_tract_df = pd.DataFrame(bike_station_to_tract.drop(columns='geometry'))
bike_station_by_tract = bike_station_to_tract_df.groupby(['GEOID'])[['id']].count()
bike_station_by_tract = bike_station_by_tract.reset_index()
bike_station_by_tract.columns = ['GEOID', 'STATIONS']

us_census_tract_df = pd.DataFrame(us_census_tract.drop(columns='geometry'))
land_area_df = us_census_tract_df[['GEOID', 'ALAND']]


acs_population = acs_population[['GEOID', 'populationE']]
acs_population['GEOID'] = acs_population['GEOID'].astype(str).str.zfill(11)
land_area_df['GEOID'] = land_area_df['GEOID'].astype(str).str.zfill(11)
bike_station_by_tract['GEOID'] = bike_station_by_tract['GEOID'].astype(str).str.zfill(11)

acs_population_with_bike = pd.merge(acs_population, land_area_df, on = 'GEOID', how = 'left')
acs_population_with_bike = pd.merge(acs_population_with_bike, bike_station_by_tract, on = 'GEOID', how = 'left')

# <codecell>

# calculate bike density
acs_population_with_bike = acs_population_with_bike.fillna(0)
density_var_area_name = 'bike_station_per_km2_' + str(analysis_year)
density_var_pop_name = 'bike_station_per_ppl_' + str(analysis_year)
acs_population_with_bike.loc[:, density_var_area_name] = 0
acs_population_with_bike.loc[:, density_var_pop_name] = 0

# Only compute density for non-zero zones
criteria = acs_population_with_bike['STATIONS'] > 0

acs_population_with_bike.loc[criteria, density_var_area_name] = \
    acs_population_with_bike.loc[criteria, 'STATIONS'] / acs_population_with_bike.loc[criteria, 'ALAND'] * (10**6) 
    
acs_population_with_bike.loc[criteria, density_var_pop_name] = \
    acs_population_with_bike.loc[criteria, 'STATIONS'] / acs_population_with_bike.loc[criteria, 'populationE']

# <codecell>

# write output
acs_population_with_bike = acs_population_with_bike[['GEOID', density_var_area_name, density_var_pop_name]]
if run_type == 2:
    acs_population_with_bike.to_csv(os.path.join('Network', output_path, 'bike_availability_' + str(analysis_year) + '.csv'), index = False)
else:
    acs_population_with_bike.to_csv(os.path.join('mode_choice_data_prep', 'input', 'bike_availability_' + str(analysis_year) + '.csv'), index = False)