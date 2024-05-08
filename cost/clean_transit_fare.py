# -*- coding: utf-8 -*-
"""
Created on Tue May  7 16:57:31 2024

@author: xiaodanxu
"""

import pandas as pd
from pandas import read_csv
import os
from os import listdir
import numpy as np

os.chdir('C:/FHWA_R2')

'''
from 0_clean_transit_costs.R
    Mode %in% c("MB", "CB") ~ "bus", #Aggregate BTS modes to Task 2 modes
    Mode  == "LR" ~ "light_rail",
    Mode %in% c("CR", "HR", "YR", "SR") ~ "commuter_rail")) %>% 
'''
transit_mode = {'bus': ["MB", "CB"],
                'rail_l': ["LR"],
                'rail_c': ["CR", "HR", "YR", "SR"]}
# load city to tract crosswalk
city_to_tract_file = 'spatial_boundary/CleanData/ZIP_COUNTY_LOOKUP_2023.csv'
city_to_tract = read_csv(city_to_tract_file)
city_to_tract = city_to_tract[['geoid', 'city', 'state']]
city_to_tract = city_to_tract.drop_duplicates(keep = 'first')

# load APTA fare data
# link: https://www.apta.com/research-technical-resources/transit-statistics/fare-database/
APTA_fare_file = 'Cost/RawData/2017-APTA-Fare-Database.xlsx'
tab_name = 'Fixed-Route Fares'
transit_fare = pd.read_excel(APTA_fare_file, sheet_name = tab_name)
transit_fare = transit_fare[['Mode Code', 'City', 'State', 'Adult Base Fare']]
transit_fare.loc[:, 'mode'] = np.nan

unique_cities = transit_fare.groupby(['City', 'State']).size().reset_index()

# assign GEMS modes
for md in transit_mode.keys():
    # print(md)
    transit_fare.loc[transit_fare['Mode Code'].isin(transit_mode[md]), 'mode'] = md

# calculate fare by city
transit_fare = transit_fare.dropna()
transit_fare_by_city = \
    transit_fare.groupby(['City', 'State', 'mode']).agg({'Adult Base Fare': ['mean', 'min', 'max']}) 
transit_fare_by_city = transit_fare_by_city.reset_index()
transit_fare_by_city.columns = ['city', 'state', 'mode', 'fare', 'mini', 'maxi']


transit_fare_by_city['city'] = transit_fare_by_city['city'].str.upper()
unique_cities_before_join = transit_fare_by_city.groupby(['city', 'state']).size().reset_index()

# assign fare to tract
transit_fare_by_tract = pd.merge(city_to_tract, transit_fare_by_city,
                                 on = ['city', 'state'], how = 'inner')

transit_fare_by_tract = \
    transit_fare_by_tract.drop_duplicates(subset = ['geoid', 'mode'], keep = 'first')
print('number of tracts with observations:')
print(len(transit_fare_by_tract.geoid.unique()))
unique_cities_after_join = transit_fare_by_tract.groupby(['city', 'state']).size().reset_index()

# format output
transit_fare_out = transit_fare_by_tract[['geoid', 'mode', 'state', 'fare', 'mini', 'maxi']]
transit_fare_out = transit_fare_out.rename(columns = {'geoid': 'tractcode'})
transit_fare_out.loc[:, 'tractcode'] = transit_fare_out.loc[:, 'tractcode'].astype(str).str.zfill(11)
transit_fare_out.loc[:, 'countycode'] = transit_fare_out.loc[:, 'tractcode'].str[0:5]

transit_fare_out.to_csv('Cost/CleanData/transit_fate_by_tract_2017.csv', index = False)

