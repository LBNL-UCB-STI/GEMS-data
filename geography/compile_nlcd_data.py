# -*- coding: utf-8 -*-
"""
Created on Tue Jan 30 09:54:07 2024

@author: xiaodanxu
"""

from pandas import read_csv
import os
from os import listdir
import pandas as pd
import matplotlib.pyplot as plt

os.chdir('C:/FHWA_R2/Land_use')
data_dir = 'CleanData'
grid_code_lookup = {'0': 'Unclassified',
                   '11': 'Water',
                    '12': 'Perennial Ice/Snow',
                    '21': 'Open Development',
                    '22': 'Low Intensity Development',
                    '23': 'Med Intensity Development',
                    '24': 'High Intensity Development',
                    '31': 'Barren Land',
                    '41': 'Deciduous Forest',
                    '42': 'Evergreen Forest',
                    '43': 'Mixed Forest',
                    '52': 'Shrub/Scrub',
                    '71': 'Grassland/Herbaceous',
                    '81': 'Pasture Land',
                    '82': 'Crop Land',
                    '90': 'Woody Wetlands',
                    '95': 'Emergent Herbaceous Wetlands'
                    }

grid_code_ak_lookup = {'0': 'Unclassified',
                   '11': 'Water',
                    '12': 'Perennial Ice/Snow',
                    '21': 'Open Development',
                    '22': 'Low Intensity Development',
                    '23': 'Med Intensity Development',
                    '24': 'High Intensity Development',
                    '31': 'Barren Land',
                    '41': 'Deciduous Forest',
                    '42': 'Evergreen Forest',
                    '43': 'Mixed Forest',
                    '51': 'Dwarf Scrub',
                    '52': 'Shrub/Scrub',
                    '71': 'Grassland/Herbaceous',
                    '72': 'Sedge/Herbaceous',
                    '74': 'Moss',
                    '81': 'Pasture Land',
                    '82': 'Crop Land',
                    '90': 'Woody Wetlands',
                    '95': 'Emergent Herbaceous Wetlands'
                    }

conus_ak_crosswalk = {'Unclassified': 'Unclassified',
                      'Water': 'Water',
                      'Perennial Ice/Snow': 'Iceland',
                      'Open Development': 'Developed Open Space',
                      'Low Intensity Development': 'Impervious Developed',
                      'Med Intensity Development': 'Impervious Developed',
                      'High Intensity Development': 'Impervious Developed',
                      'Barren Land': 'Bare Land',
                      'Deciduous Forest': 'Forest',
                      'Evergreen Forest': 'Forest',
                      'Mixed Forest': 'Forest',
                      'Shrub/Scrub': 'Shrub/Scrub',                     
                      'Grassland/Herbaceous': 'Grassland',
                      'Sedge/Herbaceous': 'Grassland',
                      'Moss': 'Grassland',
                      'Pasture Land': 'Pasture Land',
                      'Crop Land': 'Crop Land',
                      'Woody Wetlands': 'Wetland',
                      'Emergent Herbaceous Wetlands':'Wetland'
                    }

hawaii_crosswalk = {'Evergreen Forest': 'Forest', 
                    'Impervious Surface': 'Impervious Developed', 
                    'Grassland': 'Grassland', 
                    'Bare Land': 'Bare Land',	
                    'Scrub/Shrub': 'Shrub/Scrub', 
                    'Developed Open Space': 'Developed Open Space', 
                    'Open Space Developed': 'Developed Open Space',
                    'Water': 'Water',	
                    'Palustrine Emergent Wetland': 'Wetland',	
                    'Pasture/Hay': 'Pasture Land', 
                    'Cultivated': 'Crop Land',	
                    'Palustrine Forested Wetland': 'Wetland', 
                    'Estuarine Emergent Wetland':'Wetland',	
                    'Unconsolidated Shore': 'Shoreland',	
                    'Palustrine Scrub/Shrub Wetland': 'Wetland', 
                    'Palustrine Aquatic Bed': 'Shoreland',	
                    'Estuarine Forested Wetland': 'Wetland',	
                    'Estuarine Scrub/Shrub Wetland': 'Wetland', 
                    'Unclassified': 'Unclassified'
}

# <codecell>

# load and process CONUS (2020)
CONUS_NLCD = read_csv(os.path.join(data_dir, 'tract_level_land_use_no_ak.csv'))

gridcode_list = list(grid_code_lookup.keys())
print(gridcode_list)
CONUS_NLCD_long = pd.melt(CONUS_NLCD, id_vars=['GEOID'], value_vars= gridcode_list,
                          var_name = 'GRIDCODE', value_name = 'area_m2')
CONUS_NLCD_long.loc[:, 'code_name'] = CONUS_NLCD_long.loc[:, 'GRIDCODE'].map(grid_code_lookup)
CONUS_NLCD_long.loc[:, 'land_type'] = CONUS_NLCD_long.loc[:, 'code_name'].map(conus_ak_crosswalk)

CONUS_NLCD_long = CONUS_NLCD_long.groupby(['GEOID', 'land_type'])[['area_m2']].sum()
CONUS_NLCD_long = CONUS_NLCD_long.reset_index()

# <codecell>

# load and process Alaska (2010)
AK_NLCD = read_csv(os.path.join(data_dir, 'tract_level_land_use_ak.csv'))

gridcode_ak_list = list(grid_code_ak_lookup.keys())
print(gridcode_ak_list)

AK_NLCD_long = pd.melt(AK_NLCD, id_vars=['GEOID'], value_vars= gridcode_ak_list,
                          var_name = 'GRIDCODE', value_name = 'area_m2')
AK_NLCD_long.loc[:, 'code_name'] = AK_NLCD_long.loc[:, 'GRIDCODE'].map(grid_code_ak_lookup)
AK_NLCD_long.loc[:, 'land_type'] = AK_NLCD_long.loc[:, 'code_name'].map(conus_ak_crosswalk)

AK_NLCD_long = AK_NLCD_long.groupby(['GEOID', 'land_type'])[['area_m2']].sum()
AK_NLCD_long = AK_NLCD_long.reset_index()

# <codecell>

# load and process Hawaii (2010)
HI_files = ['HI_NLCD_2010_hawaii.csv', 'HI_NLCD_2010_maui.csv',
            'HI_NLCD_2010_oahu.csv', 'HI_NLCD_2010_kauai.csv',
            'HI_NLCD_2010_molokai.csv', 'HI_NLCD_2010_lanai.csv',
            'HI_NLCD_2010_niihau.csv']

HI_NLCD = pd.concat((read_csv(os.path.join(data_dir, f)) for f in HI_files), ignore_index=True)

gridcode_hi_list = list(hawaii_crosswalk.keys())
print(gridcode_hi_list)

HI_NLCD_long = pd.melt(HI_NLCD, id_vars=['GEOID'], value_vars= gridcode_hi_list,
                          var_name = 'GRIDCODE', value_name = 'area_m2')
HI_NLCD_long = HI_NLCD_long.dropna()
HI_NLCD_long.loc[:, 'land_type'] = HI_NLCD_long.loc[:, 'GRIDCODE'].map(hawaii_crosswalk)

HI_NLCD_long = HI_NLCD_long.groupby(['GEOID', 'land_type'])[['area_m2']].sum()
HI_NLCD_long = HI_NLCD_long.reset_index()


# <codecell>
US_NLCD_combined = pd.concat([CONUS_NLCD_long, AK_NLCD_long, HI_NLCD_long])
US_NLCD_combined.loc[:,  'fraction'] = US_NLCD_combined.loc[:,  'area_m2'] / \
    US_NLCD_combined.groupby('GEOID')['area_m2'].transform('sum')
US_NLCD_combined.to_csv(os.path.join(data_dir, 'processed_NLCD_data.csv'), index = False)

# <codecell>

import pygris
import geopandas as gpd

ca_tract = pygris.tracts(state = "CA", year = 2021)

ca_tract.plot()

# <codecell>
#plot ca developed land 
US_NLCD_selected = US_NLCD_combined.loc[US_NLCD_combined['land_type'] == 'Impervious Developed']
US_NLCD_selected.loc[:, 'GEOID'] = US_NLCD_selected.loc[:, 'GEOID'].astype(str).str.zfill(11)
ca_tract = ca_tract.merge(US_NLCD_selected, on ='GEOID', how = 'left')



# <codecell>


ca_tract.plot(column = 'fraction', legend=True)
plt.title('Fraction of developed land')
plt.savefig(os.path.join(data_dir, 'sample_NLCD_developed_land.png'), dpi = 300)

ca_tract.to_file(os.path.join(data_dir, "sample_ca_nlcd.geojson"), driver='GeoJSON')


# <codecell>

# spatial imputation of missing values --> use value from nearest tracts with values
select_attr = ['Impervious Developed', 'Developed Open Space']
US_NLCD_selected =\
    US_NLCD_combined.loc[US_NLCD_combined['land_type'].isin(select_attr)]
US_NLCD_wide = pd.pivot_table(US_NLCD_selected, values = 'fraction', index=['GEOID'],
                       columns = 'land_type', aggfunc="mean")
US_NLCD_wide = US_NLCD_wide.reset_index()

# <codecell>
# load census tract boundary
ct_file = '../spatial_boundary/CleanData/combined_tracts_2020.geojson'
us_census_tract = gpd.read_file(ct_file)

us_census_tract = us_census_tract.to_crs(epsg=3857)
us_census_tract_centroid = us_census_tract.centroid

us_census_tract_df = pd.DataFrame(us_census_tract.drop(columns='geometry'))
us_census_tract_centroid = pd.concat([us_census_tract_df, us_census_tract_centroid], axis = 1)
us_census_tract_centroid = gpd.GeoDataFrame(us_census_tract_centroid, geometry=0)

# <codecell>
US_NLCD_wide['GEOID'] = US_NLCD_wide['GEOID'].astype(str).str.zfill(11)
us_census_tract_centroid['GEOID'] = us_census_tract_centroid['GEOID'].astype(str).str.zfill(11)
us_census_tract_with_nlcd = us_census_tract_centroid.merge(US_NLCD_wide,
                                                          on = 'GEOID',
                                                          how = 'left')
# <codecell>
# impute from nearest location with values
# Create a spatial index

us_census_tract_no_missing = us_census_tract_with_nlcd.loc[~us_census_tract_with_nlcd['Impervious Developed'].isna()]
us_census_tract_with_missing = us_census_tract_with_nlcd.loc[us_census_tract_with_nlcd['Impervious Developed'].isna()]
sindex = us_census_tract_no_missing.sindex
us_census_tract_no_missing = us_census_tract_no_missing.reset_index()
# Find the nearest feature for each missing value
nearest_features = \
    sindex.nearest(us_census_tract_with_nlcd.geometry[us_census_tract_with_nlcd['Impervious Developed'].isna()])

# Impute the missing values with the value of the nearest feature
us_census_tract_with_missing = us_census_tract_with_missing.drop(columns = select_attr)

imputed_1 = \
    us_census_tract_no_missing['Impervious Developed'][nearest_features[1,]]
    
imputed_2 = \
        us_census_tract_no_missing['Developed Open Space'][nearest_features[1,]]
        
us_census_tract_with_missing = pd.concat([us_census_tract_with_missing.reset_index(), 
                                          imputed_1.reset_index(), 
                                          imputed_2.reset_index()], axis = 1)

us_census_tract_with_missing = us_census_tract_with_missing.drop(columns = 'index')
us_census_tract_no_missing = us_census_tract_no_missing.drop(columns = 'index')

us_census_tract_with_nlcd = pd.concat([us_census_tract_no_missing, us_census_tract_with_missing])
ca_census_tract_with_nlcd = us_census_tract_with_nlcd.loc[us_census_tract_with_nlcd['STATEFP'] == '06']
ca_census_tract_with_nlcd.to_file(os.path.join(data_dir, "sample_ca_nlcd_imputed.geojson"), driver='GeoJSON')
    
us_census_tract_with_nlcd_df = pd.DataFrame(us_census_tract_with_nlcd.drop(columns=0))
us_census_tract_with_nlcd_df = us_census_tract_with_nlcd_df[['GEOID', 'Impervious Developed', 'Developed Open Space']]

us_census_tract_with_nlcd_df.to_csv(os.path.join(data_dir, 'imputed_NLCD_data_dev_only.csv'), index = False)