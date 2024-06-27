# -*- coding: utf-8 -*-
"""
Created on Wed Jun 26 09:12:18 2024

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

run_type_desc = {1: 'mode choice', 2: 'GEMS update'}
run_type = 1
if run_type == 1:
    ct_file = 'spatial_boundary/CleanData/combined_tracts_2018.geojson'
    analysis_year = 2017
else:
    ct_file = 'spatial_boundary/CleanData/combined_tracts_2020.geojson'
    analysis_year = 2021

us_census_tract = gpd.read_file(ct_file)
data_crs = us_census_tract.crs


list_of_opportunity_files = ['AllPlacesOfWorship_4715477113142218993.geojson',
                             'Child_Care_Centers.geojson',
                             'Colleges_and_Universities_-7164156980811705283.geojson',
                             'FDIC_Insured_Banks.geojson',
                             'Hospitals.geojson',
                             'Nursing_Homes.geojson',
                             'Pharmacies_.geojson',
                             'Private_Schools_-8320420491887656674.geojson',
                             'Public_Schools_4708448687569744255.geojson',
                             'Urgent_Care_Facilities_8197288591703272850.geojson',
                             'Veterans_Health_Administration_Medical_Facilities_-2347442371999898395.geojson']

opportunity_var_name = {'AllPlacesOfWorship_4715477113142218993.geojson': 'religion',
                             'Child_Care_Centers.geojson':'childcare',
                             'Colleges_and_Universities_-7164156980811705283.geojson':'university',
                             'FDIC_Insured_Banks.geojson':'bank',
                             'Hospitals.geojson':'hospital',
                             'Nursing_Homes.geojson':'nursing_home',
                             'Pharmacies_.geojson':'pharmacy',
                             'Private_Schools_-8320420491887656674.geojson':'private_school',
                             'Public_Schools_4708448687569744255.geojson':'public_school',
                             'Urgent_Care_Facilities_8197288591703272850.geojson':'urgent_care',
                             'Veterans_Health_Administration_Medical_Facilities_-2347442371999898395.geojson':'va_medical'}


# <codecell>


def assign_points_to_tracts(points, tracts, var_name, data_crs):
    points = points.to_crs(data_crs)
    points_to_tract = points.sjoin(tracts, how="left")
    points_to_tract_df = pd.DataFrame(points_to_tract.drop(columns='geometry'))
    points_by_tract = points_to_tract_df.groupby(['GEOID'])[[var_name]].count()
    points_by_tract = points_by_tract.reset_index()
    return points_by_tract


us_census_tract_df = pd.DataFrame(us_census_tract.drop(columns='geometry'))
output_opportunities = us_census_tract_df[['GEOID', 'ALAND']]

# <codecell>

output_opportunities = us_census_tract_df[['GEOID', 'ALAND']]

for file in list_of_opportunity_files:
    
    file_path = os.path.join('Opportunity', data_path, file)
    opportunity_gdf = gpd.read_file(file_path)
    var_name = opportunity_var_name[file]
    print('generating opportunity for ' + var_name)
    opportunity_gdf.loc[:, var_name] = 1 # add a variable for counting purpose
    opportunity_by_tract = assign_points_to_tracts(opportunity_gdf, us_census_tract, var_name, data_crs)
    print(opportunity_by_tract.head(5))
    output_opportunities = pd.merge(output_opportunities, opportunity_by_tract,
                                    on = 'GEOID', how = 'left')
    
    
    # break
output_opportunities = output_opportunities.fillna(0)
# <codecell>

# load job count
if run_type == 1:
    employment_location_data = read_csv('Demand/CleanData/wac_tract_2017.csv')
else:
    employment_location_data = read_csv('Demand/CleanData/wac_tract_2021.csv')

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
employment_location_short = employment_location_short.reset_index()
employment_location_short['GEOID'] = employment_location_short['GEOID'].astype(str).str.zfill(11)


output_opportunities = pd.merge(output_opportunities, employment_location_short,
                                on = 'GEOID', how = 'left')
output_opportunities = output_opportunities.fillna(0)

# <codecell>

# save output before assigning parks
output_dir = os.path.join('Opportunity', output_path, 'opportunities_and_jobs_no_parks.csv')
output_opportunities.to_csv(output_dir, index = False)
# <codecell>

# process USA parks
us_park_file = os.path.join('Opportunity', data_path, 'USA_parks.geojson')
us_parks = gpd.read_file(us_park_file)
us_parks = us_parks.to_crs(data_crs)

# try intersection
sqmi_to_m2 = 2.59* 10**6
print(us_parks['SQMI'].sum()*sqmi_to_m2)
parks_by_tracts = gpd.overlay(us_parks, us_census_tract, how='intersection')

# <codecell>

# plot sample

sample_park_name = 'Yosemite National Park'
us_parks_sample = us_parks.loc[us_parks['NAME'] == sample_park_name]
ax = us_parks_sample.plot(column = 'OBJECTID', alpha = 0.5)

parks_by_tracts_sample = parks_by_tracts.loc[parks_by_tracts['NAME_1'] == sample_park_name]
parks_by_tracts_sample.plot(column = 'GEOID', alpha = 0.5)


parks_by_tracts.loc[:, 'park_area'] = \
    parks_by_tracts['geometry'].to_crs({'proj':'cea'}).area

print(parks_by_tracts.loc[:, 'park_area'].sum())

park_area_by_tract = \
    parks_by_tracts.groupby('GEOID')[['park_area']].sum()
park_area_by_tract = park_area_by_tract.reset_index()    
park_count_by_tract = \
        parks_by_tracts.groupby('GEOID')[['NAME_1']].count()
park_count_by_tract.columns = ['park_count']
park_count_by_tract = park_count_by_tract.reset_index()
park_summary_by_tract = pd.merge(park_area_by_tract,
                                 park_count_by_tract,
                                 on = 'GEOID', how = 'left')

# <codecell>
output_opportunities = pd.merge(output_opportunities, 
                                park_summary_by_tract,
                                on = 'GEOID', how = 'left')
output_opportunities = output_opportunities.fillna(0)
output_opportunities.loc[:, 'park_area_frac'] = \
    output_opportunities.loc[:, 'park_area'] / \
        output_opportunities.loc[:, 'ALAND']

output_dir_2 = os.path.join('Opportunity', output_path, 'opportunities_and_jobs_parks.csv')
output_opportunities.to_csv(output_dir_2, index = False)

