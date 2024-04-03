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
land_use_data = read_csv('Land_use/CleanData/imputed_NLCD_data_dev_only.csv')
urban_definition_data = read_csv('spatial_boundary/CleanData/urban_divisions_2021.csv')
path_to_emp_od = 'Demand/CleanData/OD_distance'
employment_od_data = listdir(path_to_emp_od)

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
part 2 - diversity based attributes
"""

### job residence balance--> use person as there are more non-zero values
output_demand_attributes.loc[:, 'jobs_resident_bal'] = \
    output_demand_attributes.loc[:, 'total_jobs'] / output_demand_attributes.loc[:, 'populationE']

output_demand_attributes.loc[:, 'jobs_resident_bal'] = \
    output_demand_attributes.loc[:, 'jobs_resident_bal'].fillna(0)  # NA created if 0/0, replace with 0
    
# check missing
print('total missing values in job house balance is: ')
print(output_demand_attributes.loc[:, 'jobs_resident_bal'].isnull().sum())

# check infinity
print('total infinity values in job house balance is: ')
print(sum(np.isinf(output_demand_attributes.loc[:, 'jobs_resident_bal'])))

# <codecell>
### job diversity

# calculate 8-tier employment classification
employment_location_data.loc[:, 'office_jobs'] = \
    employment_location_data.loc[:, 'naics_51'] + employment_location_data.loc[:, 'naics_52'] + \
    employment_location_data.loc[:, 'naics_53'] + employment_location_data.loc[:, 'naics_55']

employment_location_data.loc[:, 'retail_jobs'] = \
    employment_location_data.loc[:, 'naics_4445']

employment_location_data.loc[:, 'industry_jobs'] = \
    employment_location_data.loc[:, 'naics_11'] + employment_location_data.loc[:, 'naics_21'] + \
    employment_location_data.loc[:, 'naics_22'] + employment_location_data.loc[:, 'naics_23'] + \
    employment_location_data.loc[:, 'naics_3133'] + employment_location_data.loc[:, 'naics_42'] + \
    employment_location_data.loc[:, 'naics_4849'] 
    
employment_location_data.loc[:, 'service_jobs'] = \
    employment_location_data.loc[:, 'naics_54'] + employment_location_data.loc[:, 'naics_56'] + \
    employment_location_data.loc[:, 'naics_81']  
    
employment_location_data.loc[:, 'recreation_jobs'] = \
    employment_location_data.loc[:, 'naics_71'] + employment_location_data.loc[:, 'naics_72']

employment_location_data.loc[:, 'education_jobs'] = \
    employment_location_data.loc[:, 'naics_61']
    
employment_location_data.loc[:, 'healthcare_jobs'] = \
    employment_location_data.loc[:, 'naics_62']

employment_location_data.loc[:, 'government_jobs'] = \
    employment_location_data.loc[:, 'naics_92']

tier_list = ['office_jobs', 'retail_jobs', 'industry_jobs', 'service_jobs',
             'recreation_jobs',  'education_jobs', 'healthcare_jobs', 'government_jobs']

employment_location_data.loc[:, 'total_jobs'] = \
    employment_location_data.loc[:, 'total_jobs'].fillna(0)
    
employment_location_data.loc[:, tier_list] = \
    employment_location_data[tier_list].div(employment_location_data['total_jobs'], axis=0)
    
entropy_list = []
for tier in tier_list:
    # print('calculate entropy for ' + tier)
    e_var = 'e_' + tier
    entropy_list.append(e_var)
    employment_location_data.loc[:, e_var] = 0
    
    # only calculate entropy for non-zero values
    non_zero_id = (employment_location_data[tier] > 0)
    employment_location_data.loc[non_zero_id, e_var] = \
        -employment_location_data.loc[non_zero_id, tier] * \
        np.log(employment_location_data.loc[non_zero_id, tier])

# if the tract has 0 employment, it will not have job diversity as outcome -> diversity = N/A
employment_location_data.loc[:, 'job_diversity'] = np.nan

# if diversity = 0, it means the tract only has 1 industry
non_zero_zone = (employment_location_data['total_jobs'] > 0)
employment_location_data.loc[non_zero_zone, 'job_diversity'] = \
        employment_location_data.loc[non_zero_zone, entropy_list].sum(axis = 1) / \
        np.log(8) 

employment_diversity = employment_location_data[['GEOID', 'job_diversity']]        
output_demand_attributes = pd.merge(output_demand_attributes,
                                    employment_diversity,
                                    on = 'GEOID', how = 'left')

# check missing
print('total missing values in job diversity is: ')
print(output_demand_attributes.loc[:, 'job_diversity'].isnull().sum())

# check infinity
print('total infinity values in job diversity is: ')
print(sum(np.isinf(output_demand_attributes.loc[:, 'job_diversity'])))


# <codecell>

### job sink magnitude and trip distance distribution

## update April 01, 2024 --split bin 4 at 50 mile per Mona's request
od_attribute_exist = 0 # if 0, execute the data generation, if 1, load existing output

out_od_attributes = None # create empty data frame to hold generated attributes

dist_bin = [-1, 1.3, 3, 8, 50, 150, 5000] 
dist_bin_label = ['jobs 0-1.3 miles', 'jobs 1.3-3 miles', 'jobs 3-8 miles', 
                  'jobs 8-50 miles', 'jobs 50-150 miles', 'remote jobs']

if od_attribute_exist == 0: # run heavy comutation to generate od attributes
    od_dist_by_tract = None  # create empty data frame to hold input data
    trip_by_home_tract = None  # create empty data frame to hold input data
    trip_by_work_tract = None  # create empty data frame to hold input data
    for data in employment_od_data:
        print('processing od data ' + data)
        od_data = read_csv(os.path.join(path_to_emp_od, data))
        
        # trip distance distribution
        od_data.loc[:, 'dist_bin'] = \
            pd.cut(od_data.loc[:, 'distance'], bins = dist_bin,
           labels = dist_bin_label)
        od_data_by_dist = pd.pivot_table(od_data, 
                                           values = 'S000',
                                           index = 'w_tract', 
                                           columns = 'dist_bin',
                                           aggfunc = "sum")
        od_data_by_dist = od_data_by_dist.reset_index()
        od_data_by_dist = od_data_by_dist.fillna(0)
        od_data_by_dist.loc[:, 'total_emp'] = \
            od_data_by_dist.loc[:, dist_bin_label].sum(axis = 1)
        od_data_by_dist = od_data_by_dist.loc[od_data_by_dist['total_emp'] > 0 ] 
        # fraction only available to zones with non-zero jobs
        
        # calculate fraction of trips by distance bin
        od_data_by_dist.loc[:, dist_bin_label] = od_data_by_dist.loc[:, dist_bin_label].div(
            od_data_by_dist.loc[:, 'total_emp'], axis = 0)
        od_dist_by_tract = pd.concat([od_dist_by_tract, od_data_by_dist])
        
        
        # trip sink magnitude (aggregation performed after concat to account for work/home in different states)
        od_by_home = od_data.groupby('h_tract')[['S000']].sum()
        
        od_by_home = od_by_home.reset_index()
        od_by_home.columns = ['GEOID', 'jobs_by_home']
        trip_by_home_tract = pd.concat([trip_by_home_tract, od_by_home])
        
        od_by_work = od_data.groupby('w_tract')[['S000']].sum()
        od_by_work = od_by_work.reset_index()
        od_by_work.columns = ['GEOID', 'jobs_by_work']
        trip_by_work_tract = pd.concat([trip_by_work_tract, od_by_work])
        # break
    
    print('if O-D data contains duplicated geoid:')
    print(od_dist_by_tract['w_tract'].duplicated().any())
    trip_by_home_tract = trip_by_home_tract.groupby('GEOID')[['jobs_by_home']].sum()
    trip_by_home_tract = trip_by_home_tract.reset_index()
    
    trip_by_work_tract = trip_by_work_tract.groupby('GEOID')[['jobs_by_work']].sum()
    trip_by_work_tract = trip_by_work_tract.reset_index()
    
    out_od_attributes = od_dist_by_tract
    out_od_attributes = out_od_attributes.rename(columns = {'w_tract': 'GEOID'})
    
    out_od_attributes = pd.merge(out_od_attributes, trip_by_home_tract, 
                                 on = 'GEOID', how = 'left')
    
    out_od_attributes = pd.merge(out_od_attributes, trip_by_work_tract, 
                                 on = 'GEOID', how = 'left')
    
    out_od_attributes = out_od_attributes.fillna(0)
    
    out_od_attributes.loc[:, 'job_sink_mag'] = \
        out_od_attributes.loc[:, 'jobs_by_work'] / out_od_attributes.loc[:, 'jobs_by_home']
    
    out_od_attributes.to_csv('Demand/CleanData/lehd_od_trip_characteristics_v2.csv', index = False)
else: # load pre-generated data
    out_od_attributes = read_csv('Demand/CleanData/lehd_od_trip_characteristics_v2.csv')

# <codecell>
# append OD characteristics to output metrics
var_list = ['GEOID', 'jobs 0-1.3 miles', 'jobs 1.3-3 miles', 'jobs 3-8 miles',
            'jobs 8-50 miles', 'jobs 50-150 miles', 'remote jobs', 'job_sink_mag']
out_od_attributes_short = out_od_attributes[var_list]

output_demand_attributes = pd.merge(output_demand_attributes,
                                    out_od_attributes_short,
                                    on = 'GEOID', how = 'left')

# check missing

for var in var_list:
    if var == 'GEOID':
        continue
    else:
        print('total missing values in ' + var)
        print(output_demand_attributes.loc[:, var].isnull().sum())
        
        # check infinity
        print('total infinity values in ' + var)
        print(sum(np.isinf(output_demand_attributes.loc[:, var])))
        
# <codecell>

### land use and land cover attributes
var_list = ['Impervious Developed', 'Developed Open Space']
# land_use_data_selected = land_use_data.loc[land_use_data['land_type'].isin(var_list)]
# land_use_data_wide = pd.pivot_table(land_use_data_selected,
#                                     values = 'fraction', index = ['GEOID'], 
#                                     columns = 'land_type', 
#                                     aggfunc = "sum")
# land_use_data_wide = land_use_data_wide.fillna(0)
# land_use_data_wide = land_use_data_wide.reset_index()
# land_use_data_wide.loc[:, 'Agriculture Land'] = land_use_data_wide.loc[:, 'Pasture Land'] + \
#     land_use_data_wide.loc[:, 'Crop Land']

# land_use_data_wide = \
#     land_use_data_wide[['GEOID', 'Impervious Developed', 'Developed Open Space']]

#land_use_data_wide
output_demand_attributes = pd.merge(output_demand_attributes,
                                    land_use_data,
                                    on = 'GEOID', how = 'left')

# check missing

for var in var_list:
    if var == 'GEOID':
        continue
    else:
        print('total missing values in ' + var)
        print(output_demand_attributes.loc[:, var].isnull().sum())
        
        # check infinity
        print('total infinity values in ' + var)
        print(sum(np.isinf(output_demand_attributes.loc[:, var])))
        
# <codecell>

 ### land use and land cover attributes       
var_list = ['GEOID', 'census_urban_area']
urban_definition_data_select = urban_definition_data[var_list]

output_demand_attributes = pd.merge(output_demand_attributes,
                                    urban_definition_data_select,
                                    on = 'GEOID', how = 'left')

# check missing

for var in var_list:
    if var == 'GEOID':
        continue
    else:
        print('total missing values in ' + var)
        print(output_demand_attributes.loc[:, var].isnull().sum())
        
        # check infinity
        print('total infinity values in ' + var)
        print(sum(np.isinf(output_demand_attributes.loc[:, var])))
# <codecell>

# Log April 01, 2024 - generate V2 of the data after biennual review
output_demand_attributes.to_csv('Demand/CleanData/microtype_inputs_demand_V2.csv', index = False)

