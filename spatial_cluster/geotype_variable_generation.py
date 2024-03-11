# -*- coding: utf-8 -*-
"""
Created on Mon Mar 11 10:42:49 2024

@author: xiaodanxu
"""

# load environment
import pandas as pd
from pandas import read_csv
import os
from os import listdir
import numpy as np
import matplotlib.pyplot as plt
import geopandas as gpd
from shapely.geometry import box

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
geotype_boundary_df = read_csv('spatial_boundary/CleanData/combined_geotype_unit_' + census_year + '.csv')
spatial_crosswalk = read_csv('spatial_boundary/CleanData/cleaned_lodes8_crosswalk_with_ID.csv')
employment_location_data = read_csv('Demand/CleanData/wac_tract_2021.csv')
land_use_data = read_csv('Land_use/CleanData/processed_NLCD_data.csv')
path_to_emp_od = 'Demand/CleanData/OD_distance'
employment_od_data = listdir(path_to_emp_od)

# load spatial boundary
gt_file = 'spatial_boundary/CleanData/combined_geotype_unit_' + census_year + '.geojson'
us_geotype_boudary = gpd.read_file(gt_file)

# <codecell>

# chop the map and only include CONUS
polygon = box(-130, 10, -66.96466, 50.365162)
us_geotype_boudary_clipped = us_geotype_boudary.clip(polygon)
us_geotype_boudary_clipped.loc[:,'spatial_id'] = \
us_geotype_boudary_clipped.loc[:,'spatial_id'].astype(int)
us_geotype_boudary_clipped.plot()

# <codecell>

spatial_crosswalk = spatial_crosswalk.rename(columns = {'trct': 'GEOID'})    

#create output dataframe
output_geotype_attributes = geotype_boundary_df[['spatial_id', 'is_cbsa', 'ALAND', 'AWATER', 'pct_water']]
output_geotype_attributes.loc[:,'spatial_id'] = output_geotype_attributes.loc[:,'spatial_id'].astype(int)
# remove tracts with all waters
output_geotype_attributes = \
    output_geotype_attributes.loc[output_geotype_attributes['pct_water'] < 1]
    
tract_attributes = boundary_and_area[['GEOID', 'ALAND', 'pct_water']]

# remove tracts with all waters
tract_attributes = \
    tract_attributes.loc[tract_attributes['pct_water'] < 1]

# <codecell>
# compute fraction of agriculture land
agri_land_type = ['Crop Land', 'Pasture Land']
land_use_data_agriculture = pd.pivot_table(land_use_data, 
                                           values='area_m2', index=['GEOID'],
                                           columns=['land_type'], aggfunc="sum")

land_use_data_agriculture = land_use_data_agriculture.fillna(0)
land_use_data_agriculture.loc[:, 'total_land'] = land_use_data_agriculture.sum(axis = 1)
land_use_data_agriculture.loc[:, 'agri_land'] = land_use_data_agriculture.loc[:, agri_land_type].sum(axis = 1)

land_use_data_agriculture = land_use_data_agriculture.reset_index()
land_use_data_agriculture = land_use_data_agriculture[['GEOID', 'total_land', 'agri_land']]


spatial_crosswalk_with_agri = pd.merge(spatial_crosswalk, land_use_data_agriculture,
                                       on = 'GEOID', how = 'inner')
agriculture_by_spatial_id = \
    spatial_crosswalk_with_agri.groupby(['spatial_id'])[['total_land', 'agri_land']].sum()
agriculture_by_spatial_id.loc[:, 'agriculture_frac'] = \
    agriculture_by_spatial_id.loc[:, 'agri_land'] / agriculture_by_spatial_id.loc[:, 'total_land']

agriculture_by_spatial_id = agriculture_by_spatial_id.reset_index()

agriculture_by_spatial_id = agriculture_by_spatial_id[['spatial_id', 'agriculture_frac']]


agriculture_by_spatial_id.loc[:,'spatial_id'] = agriculture_by_spatial_id.loc[:,'spatial_id'].astype(int)

us_geotype_boudary_to_plot = us_geotype_boudary_clipped.merge(agriculture_by_spatial_id,
                                                              on = 'spatial_id', how = 'left')

us_geotype_boudary_to_plot.plot(column = 'agriculture_frac', legend=True, 
                                cmap= 'viridis',
                                legend_kwds = {'shrink': 0.6})
plt.title('Fraction of agriculture land')
plt.savefig('Demand/Figures/geotype_fraction_of_agri_land.png', dpi= 300, bbox_inches = 'tight')


output_geotype_attributes = pd.merge(output_geotype_attributes, agriculture_by_spatial_id,
                                     on = 'spatial_id', how = 'left')

# <codecell>

# generate regional industry mix

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

out_var_list = tier_list
out_var_list.append('total_jobs')

employment_location_data = employment_location_data.set_index('GEOID')
employment_location_short = employment_location_data[out_var_list]

# <codecell>
spatial_crosswalk_with_industry = pd.merge(spatial_crosswalk, employment_location_short,
                                       left_on = 'GEOID', right_index = True, how = 'inner')
industry_by_spatial_id = \
    spatial_crosswalk_with_industry.groupby(['spatial_id'])[out_var_list].sum()

industry_by_spatial_id = industry_by_spatial_id.reset_index()    
tier_list = ['office_jobs', 'retail_jobs', 'industry_jobs', 'service_jobs',
             'recreation_jobs',  'education_jobs', 'healthcare_jobs', 'government_jobs']

for attr in tier_list:
    industry_by_spatial_id.loc[:, attr] = \
        industry_by_spatial_id.loc[:, attr] / industry_by_spatial_id.loc[:, 'total_jobs']

industry_by_spatial_id['spatial_id'] = industry_by_spatial_id['spatial_id'].astype(int)
us_geotype_boudary_to_plot = us_geotype_boudary_clipped.merge(industry_by_spatial_id,
                                                              on = 'spatial_id', how = 'left')
for attr in tier_list:
    us_geotype_boudary_to_plot.plot(column = attr, 
                                    legend=True, 
                                    cmap= 'viridis',
                                    legend_kwds = {'shrink': 0.6})
    plt.title('Fraction of ' + attr)
    plt.savefig('Demand/Figures/geotype_fraction_of_' + attr + '.png', dpi= 300, bbox_inches = 'tight')
    plt.show()


import seaborn as sns
corr_mat = np.round(industry_by_spatial_id[tier_list].corr(),2)
sns.heatmap(corr_mat, vmin=-1, vmax=1, annot = True, cmap = 'coolwarm')
plt.xticks(rotation =60)
plt.savefig('Demand/Figures/industry_correlation.png', dpi= 300, bbox_inches = 'tight')
plt.show()

# <codecell>
output_geotype_attributes = pd.merge(output_geotype_attributes, industry_by_spatial_id,
                                      on = 'spatial_id', how = 'left')

# <codecell>
output_geotype_attributes.loc[:, 'job_density'] = \
    output_geotype_attributes.loc[:, 'total_jobs']/(unit_converter['m2_to_acre'] * \
                                                    output_geotype_attributes.loc[:, 'ALAND'])