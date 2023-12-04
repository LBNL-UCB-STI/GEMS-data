# -*- coding: utf-8 -*-
"""
Created on Mon Dec  4 10:50:11 2023

Module: Generate demand-related input attributes for new micro-geotype
@author: xiaodanxu
"""

# load environment
import pandas as pd
from pandas import read_csv
import os
from os import listdir
import numpy as np

# define constant
unit_converter = {
    # distance
    'meter_to_mile': 0.000621371,
    # area
    'm2_to_acre':0.000247105,
    'acre_to_sqmile': 0.0015625,   
    }


# define data path
os.chdir('C:/FHWA_R2')


# load inputs -- the data needs to adopt the same geographic resolution for census tracts
census_year = '2020'
boundary_and_area = read_csv('spatial_boundary/CleanData/combined_tracts_' + census_year + '.csv')
population_data = read_csv('Demography/CleanData/acs_data_tracts_112023.csv')
employment_location_data = read_csv('Demand/CleanData/wac_tract_2021.csv')
employment_od_data = listdir('Demand/CleanData/OD_distance')

#create output dataframe
output_demand_attributes = boundary_and_area[['GEOID', 'ALAND', 'pct_water']]

# remove tracts with all waters
output_demand_attributes = \
    output_demand_attributes.loc[output_demand_attributes['pct_water'] < 1]
    
# <codecell>
"""
part 1 - density based attributes
"""
output_demand_attributes.loc[:, 'land_area_acre'] = \
    output_demand_attributes.loc[:, 'ALAND'] * unit_converter['m2_to_acre']

### pop density
population_data_short = population_data.loc[:, ['GEOID', 'populationE', 'householdsE']] # In variable name, E means estimate, M means MOE
output_demand_attributes = pd.merge(output_demand_attributes,
                                    population_data_short,
                                    on = 'GEOID', how = 'left')

output_demand_attributes.loc[:, 'pop_per_acre'] = \
    output_demand_attributes.loc[:, 'populationE'] / output_demand_attributes.loc[:, 'land_area_acre']
    
# check missing
print('total missing values in pop density is: ')
print(output_demand_attributes.loc[:, 'pop_per_acre'].isnull().sum())

# check infinity
print('total infinity values in pop density is: ')
print(sum(np.isinf(output_demand_attributes.loc[:, 'pop_per_acre'])))



### job density
employment_location_short = employment_location_data[['GEOID', 'total_jobs']]

output_demand_attributes = pd.merge(output_demand_attributes,
                                    employment_location_short,
                                    on = 'GEOID', how = 'left')

# fill missing with 0
output_demand_attributes.loc[:, 'total_jobs'] = \
    output_demand_attributes.loc[:, 'total_jobs'].fillna(0)   
    
output_demand_attributes.loc[:, 'jobs_per_acre'] = \
    output_demand_attributes.loc[:, 'total_jobs'] / output_demand_attributes.loc[:, 'land_area_acre']
    
# check missing
print('total missing values in job density is: ')
print(output_demand_attributes.loc[:, 'jobs_per_acre'].isnull().sum())

# check infinity
print('total infinity values in job density is: ')
print(sum(np.isinf(output_demand_attributes.loc[:, 'jobs_per_acre'])))

# <codecell>
"""
part 1 - diversity based attributes
"""

# job housing balance--> use person as there more non-zero values
output_demand_attributes.loc[:, 'jobs_house_bal'] = \
    output_demand_attributes.loc[:, 'total_jobs'] / output_demand_attributes.loc[:, 'populationE']

output_demand_attributes.loc[:, 'jobs_house_bal'] = \
    output_demand_attributes.loc[:, 'jobs_house_bal'].fillna(0)  # NA created if 0/0, replace with 0
# check missing
print('total missing values in job house balance is: ')
print(output_demand_attributes.loc[:, 'jobs_house_bal'].isnull().sum())

# check infinity
print('total infinity values in job house balance is: ')
print(sum(np.isinf(output_demand_attributes.loc[:, 'jobs_house_bal'])))