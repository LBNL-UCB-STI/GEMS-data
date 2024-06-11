# -*- coding: utf-8 -*-
"""
Created on Fri May 17 09:37:09 2024

@author: xiaodanxu
"""

import pandas as pd
from pandas import read_csv
import os
from os import listdir
import numpy as np

os.chdir('C:/FHWA_R2')

city_to_tract_file = 'spatial_boundary/CleanData/ZIP_COUNTY_LOOKUP_2023.csv'
city_to_tract = read_csv(city_to_tract_file)
city_to_tract = city_to_tract[['geoid', 'city', 'state']]
city_to_tract = city_to_tract.drop_duplicates(keep = 'first')


# processing parking fee
parking_fee_file = 'Cost/RawData/parkopedia_cleaned_name.csv'
parking_fee = read_csv(parking_fee_file)

parking_fee['city'] = parking_fee['city'].str.upper()
city_list_before = parking_fee['city'].unique()

parking_fee_by_tract = pd.merge(city_to_tract, parking_fee,
                                 on = ['city', 'state'], how = 'inner')
parking_fee_by_tract = \
    parking_fee_by_tract.drop_duplicates(subset = 'geoid', keep = 'first')

parking_fee_by_tract['parking'] /= 2 # convert to per hour cost     
city_list_after = parking_fee_by_tract['city'].unique()

diff_cities = set(city_list_before) - set(city_list_after)

parking_fee_by_tract = parking_fee_by_tract.rename(columns = {'geoid': 'tractcode'})
parking_fee_by_tract.loc[:, 'tractcode'] = parking_fee_by_tract.loc[:, 'tractcode'].astype(str).str.zfill(11)
parking_fee_by_tract.to_csv('Cost/CleanData/parking_tract_2017.csv', index = False)

# <codecell>
# processing Uber fee
uber_fee_file = 'Cost/RawData/Uber_fare_cleaned_name.csv'
uber_fee = read_csv(uber_fee_file)

uber_fee = uber_fee.rename(columns = {'City': 'city', 'State': 'state'})
uber_fee['city'] = uber_fee['city'].str.upper()
city_list_before = uber_fee['city'].unique()

uber_fee_by_tract = pd.merge(city_to_tract, uber_fee,
                                 on = ['city', 'state'], how = 'inner')
uber_fee_by_tract = uber_fee_by_tract.drop_duplicates(subset = 'geoid', keep = 'first')
city_list_after = uber_fee_by_tract['city'].unique()

diff_cities = set(city_list_before) - set(city_list_after)

uber_fee_by_tract = uber_fee_by_tract.rename(columns = {'geoid': 'tractcode'})
uber_fee_by_tract.loc[:, 'tractcode'] = uber_fee_by_tract.loc[:, 'tractcode'].astype(str).str.zfill(11)

uber_fee_by_tract.to_csv('Cost/CleanData/uber_fare_tract_2017.csv', index = False)