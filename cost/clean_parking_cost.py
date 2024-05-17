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

parking_fee_file = 'Cost/RawData/parkopedia_cleaned_name.csv'
parking_fee = read_csv(parking_fee_file)

parking_fee['city'] = parking_fee['city'].str.upper()
city_list_before = parking_fee['city'].unique()

parking_fee_by_tract = pd.merge(city_to_tract, parking_fee,
                                 on = ['city', 'state'], how = 'inner')
city_list_after = parking_fee_by_tract['city'].unique()

diff_cities = set(city_list_before) - set(city_list_after)

parking_fee_by_tract.to_csv('Cost/CleanData/parking_tract_2017.csv', index = False)