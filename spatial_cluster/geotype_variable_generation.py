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
network_data = read_csv('Network/CleanData/network_microtype_metrics.csv')
microtype_results_data = read_csv('Demand/Results/clustering_outputs_with_raw_data.csv')

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
output_geotype_attributes['spatial_id'] = output_geotype_attributes['spatial_id'].astype(int)
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
plt.xticks(rotation = 60)
plt.savefig('Demand/Figures/industry_correlation.png', dpi= 300, bbox_inches = 'tight')
plt.show()

# <codecell>
output_geotype_attributes = pd.merge(output_geotype_attributes, industry_by_spatial_id,
                                      on = 'spatial_id', how = 'left')

output_geotype_attributes.loc[:, 'job_density'] = \
    output_geotype_attributes.loc[:, 'total_jobs']/(unit_converter['m2_to_acre'] * \
                                                    output_geotype_attributes.loc[:, 'ALAND'])
        
# <codecell>
### commute dispersion (HHI)

od_attribute_exist = 0 # if 0, execute the data generation, if 1, load existing output

out_od_attributes = None # create empty data frame to hold generated attributes
network_data_short = network_data[['tract', 'laneMiles']]
if od_attribute_exist == 0: # run heavy comutation to generate od attributes
    od_dist_by_tract = None  # create empty data frame to hold input data
    trip_by_home_tract = None  # create empty data frame to hold input data
    trip_by_work_tract = None  # create empty data frame to hold input data
    for data in employment_od_data:
        print('processing od data ' + data)
        od_data = read_csv(os.path.join(path_to_emp_od, data))
        od_data = od_data.loc[od_data['distance'] > 150 ] # exclude long-dist commute
        od_data.loc[:, 'PMT'] = od_data.loc[:, 'S000'] * od_data.loc[:, 'distance']
        od_data = od_data.groupby(['w_tract'])[['PMT']].sum()
        od_data = od_data.reset_index()
        out_od_attributes = pd.concat([out_od_attributes, od_data])

        # break
# in case w_tract appears in multiple files, group output again
out_od_attributes = out_od_attributes.groupby(['w_tract'])[['PMT']].sum() 
out_od_attributes = out_od_attributes.reset_index() 

# <codecell>
 
# calculate SI value
out_od_attributes = out_od_attributes.rename(columns = {'w_tract': 'GEOID'})
network_data_short = network_data_short.rename(columns = {'tract': 'GEOID'})
out_od_with_lm = pd.merge(out_od_attributes, network_data_short, 
                          on = 'GEOID', how = 'left')
out_od_with_lm = out_od_with_lm.dropna()
out_od_with_lm.loc[:, 'pmt_per_lm'] = out_od_with_lm.loc[:, 'PMT'] / \
    out_od_with_lm.loc[:, 'laneMiles']

spatial_crosswalk_with_hhi = pd.merge(spatial_crosswalk, out_od_with_lm,
                                       on = 'GEOID', how = 'inner')
spatial_crosswalk_with_hhi['pmt_per_lm'].fillna(0, inplace = True)
spatial_crosswalk_with_hhi.loc[:, 'SI'] = \
    spatial_crosswalk_with_hhi.loc[:, 'pmt_per_lm'] / \
        spatial_crosswalk_with_hhi.groupby('spatial_id')['pmt_per_lm'].transform('sum')
spatial_crosswalk_with_hhi.loc[:, 'HHI'] = np.square(spatial_crosswalk_with_hhi.loc[:, 'SI'])
hhi_by_spatial_id = \
    spatial_crosswalk_with_hhi.groupby(['spatial_id'])[['HHI', 'laneMiles']].sum()


hhi_by_spatial_id = hhi_by_spatial_id.reset_index()  
hhi_by_spatial_id.loc[:,'spatial_id'] = hhi_by_spatial_id.loc[:,'spatial_id'].astype(int)


us_geotype_boudary_to_plot = \
    us_geotype_boudary_clipped.merge(hhi_by_spatial_id,
                                     on = 'spatial_id', how = 'left')

us_geotype_boudary_to_plot.plot(column = 'HHI', legend=True, 
                                cmap= 'viridis',
                                legend_kwds = {'shrink': 0.6})
plt.title('Commute dispersion index')
plt.savefig('Demand/Figures/geotype_dispersion_index.png', dpi= 300, bbox_inches = 'tight')
hhi_by_spatial_id.to_csv('Demand/CleanData/geotype_hhi_lm.csv', index = False)

# <codecell>
output_geotype_attributes = pd.merge(output_geotype_attributes, hhi_by_spatial_id,
                                      on = 'spatial_id', how = 'left')

output_geotype_attributes.loc[:, 'lm_density'] = \
    output_geotype_attributes.loc[:, 'laneMiles']/(unit_converter['m2_to_acre'] * \
                                                    output_geotype_attributes.loc[:, 'ALAND'])

# <codecell>
lm_to_plot = output_geotype_attributes[['spatial_id', 'lm_density', 'job_density']]
us_geotype_boudary_to_plot = \
    us_geotype_boudary_clipped.merge(lm_to_plot,
                                     on = 'spatial_id', how = 'left')
us_geotype_boudary_to_plot.plot(column = 'lm_density', legend=True, 
                                cmap= 'viridis',
                                legend_kwds = {'shrink': 0.6})
plt.title('Lane mile per acre')
plt.savefig('Demand/Figures/geotype_lm_per_acre.png', dpi= 300, bbox_inches = 'tight')

us_geotype_boudary_to_plot.plot(column = 'job_density', legend=True, 
                                cmap= 'viridis',
                                legend_kwds = {'shrink': 0.6})
plt.title('Jobs per acre')
plt.savefig('Demand/Figures/geotype_job_per_acre.png', dpi= 300, bbox_inches = 'tight')


# <codecell>

# fraction of microtypes
microtype_results_selected = \
microtype_results_data[['GEOID', 'ALAND', 'populationE', 'demand_microtype_comb']]
microtype_results_selected = microtype_results_selected.rename(columns = {'populationE': 'total_person'})
spatial_crosswalk_with_microtype = pd.merge(spatial_crosswalk, 
                                            microtype_results_selected,
                                       on = 'GEOID', how = 'inner')
unique_types = microtype_results_selected.demand_microtype_comb.unique()

microtype_by_spatial_id = pd.pivot_table(spatial_crosswalk_with_microtype,
                                         values = 'ALAND', index = 'spatial_id',
                                         columns = 'demand_microtype_comb', aggfunc="sum")
microtype_by_spatial_id.fillna(0, inplace = True)
microtype_by_spatial_id.loc[:, 'total'] = microtype_by_spatial_id.sum(axis = 1)

for mt in unique_types:
    microtype_by_spatial_id.loc[:, mt] /= microtype_by_spatial_id.loc[:, 'total']
microtype_by_spatial_id = microtype_by_spatial_id.reset_index()
microtype_by_spatial_id = microtype_by_spatial_id.drop(columns = 'total')

pop_density_by_spatial_id = \
    spatial_crosswalk_with_microtype.groupby('spatial_id')[['total_person', 'ALAND']].sum()
pop_density_by_spatial_id = pop_density_by_spatial_id.reset_index()

pop_density_by_spatial_id.loc[:, 'pop_density'] = \
    pop_density_by_spatial_id.loc[:, 'total_person']/(unit_converter['m2_to_acre'] * \
                                                    pop_density_by_spatial_id.loc[:, 'ALAND'])
pop_density_by_spatial_id = pop_density_by_spatial_id[['spatial_id', 'pop_density']]
microtype_by_spatial_id = pd.merge(microtype_by_spatial_id, pop_density_by_spatial_id,
                                   on = 'spatial_id', how = 'left')

# <codecell>
us_geotype_boudary_to_plot = \
    us_geotype_boudary_clipped.merge(microtype_by_spatial_id,
                                     on = 'spatial_id', how = 'left')
for mt in unique_types:
    us_geotype_boudary_to_plot.plot(column = mt, 
                                    legend=True, 
                                    cmap= 'viridis',
                                    legend_kwds = {'shrink': 0.6})
    plt.title('Fraction of ' + mt)
    plt.savefig('Demand/Figures/geotype_fraction_of_microtype_' + mt + '.png', dpi= 300, bbox_inches = 'tight')
    plt.show()    

us_geotype_boudary_to_plot.plot(column = 'pop_density', legend=True, 
                                cmap= 'viridis',
                                legend_kwds = {'shrink': 0.6})
plt.title('Population per acre')
plt.savefig('Demand/Figures/geotype_pop_per_acre.png', dpi= 300, bbox_inches = 'tight')
# <codecell>
output_geotype_attributes = pd.merge(output_geotype_attributes, microtype_by_spatial_id,
                                      on = 'spatial_id', how = 'left')
# <codecell>

# writing output
output_geotype_attributes[['HHI']] = \
    output_geotype_attributes[['HHI']].fillna(1)
output_geotype_attributes[['lm_density']] = \
    output_geotype_attributes[['lm_density']].fillna(method = 'ffill')
output_geotype_attributes.to_csv('Demand/CleanData/geotype_inputs.csv', index = False)
