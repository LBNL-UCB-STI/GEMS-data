# -*- coding: utf-8 -*-
"""
Created on Fri Jun 14 09:52:06 2024

@author: xiaodanxu
"""

# load environment
import pandas as pd
from pandas import read_csv
import os
from os import listdir
import numpy as np

# define data path
os.chdir('C:/FHWA_R2')

path_to_typology = 'Network/CleanData/'
path_to_demand = 'Demand/Results/'
path_to_crosswalk = 'spatial_boundary/CleanData/'



census_tract_crosswalk_file = 'census_tract_crosswalk_2010_2020.csv'
spatial_id_crosswalk_file = 'cleaned_lodes8_crosswalk_with_ID.csv'
spatial_id_crosswalk_2010_file = 'us_xwalk_tract_2017_withID.csv'

typology_unimputed_2020_file = 'microtype_geotype_output_2020_noimp_V1.csv'
typology_unimputed_2010_file = 'microtype_geotype_output_2010_noimp_V1.csv'
#partition_results_file = 'final_partition_results.csv'

census_tract_crosswalk = read_csv(os.path.join(path_to_crosswalk, census_tract_crosswalk_file))
spatial_id_crosswalk = read_csv(os.path.join(path_to_crosswalk, spatial_id_crosswalk_file))
spatial_id_crosswalk_2010 = read_csv(os.path.join(path_to_crosswalk, spatial_id_crosswalk_2010_file))

typology_unimputed_2020 = read_csv(os.path.join(path_to_typology, typology_unimputed_2020_file))
typology_unimputed_2010 = read_csv(os.path.join(path_to_typology, typology_unimputed_2010_file))


#partition_results = read_csv(os.path.join(path_to_typology, partition_results_file))
# spatial_id_crosswalk = spatial_id_crosswalk[['spatial_id', 'trct']]

# geotype_combined = pd.merge(spatial_id_crosswalk, geotype_combined, 
#                             on = 'spatial_id', how = 'left')

# <codecell>
typology_unimputed_2020 = typology_unimputed_2020.drop(columns = ['spatial_id'])
typology_unimputed_2020['GEOID'] = typology_unimputed_2020['GEOID'].astype(str).str.zfill(11)

spatial_id_crosswalk_short = spatial_id_crosswalk[['spatial_id', 'trct']]
spatial_id_crosswalk_short = spatial_id_crosswalk_short.rename(columns = {'trct': 'GEOID'})

spatial_id_crosswalk_short['GEOID'] = spatial_id_crosswalk_short['GEOID'].astype(str).str.zfill(11)

typology_unimputed_2020 = pd.merge(spatial_id_crosswalk_short, typology_unimputed_2020, 
                            on = 'GEOID', how = 'left')

typology_unimputed_2020 = typology_unimputed_2020.sort_values('GEOID', ascending = True)
print('total missing before imputation')
print(typology_unimputed_2020.isna().sum())

#impute rural and urban respectively
micro_geotype_combined_imputed = typology_unimputed_2020.copy()
micro_geotype_combined_imputed.loc[:, 'geotype'] = \
    micro_geotype_combined_imputed.loc[:, 'geotype'].fillna(method = 'ffill')

micro_geotype_urban_imputed = \
    micro_geotype_combined_imputed.loc[micro_geotype_combined_imputed['geotype'].isin(['CBSA_1', 'CBSA_2'])]
micro_geotype_urban_imputed = micro_geotype_urban_imputed.fillna(method = 'ffill')
# micro_geotype_urban_imputed = micro_geotype_urban_imputed.fillna(method = 'bfill')
micro_geotype_rural_imputed = \
    micro_geotype_combined_imputed.loc[~micro_geotype_combined_imputed['geotype'].isin(['CBSA_1', 'CBSA_2'])]
micro_geotype_rural_imputed = micro_geotype_rural_imputed.fillna(method = 'ffill')
# micro_geotype_rural_imputed = micro_geotype_rural_imputed.fillna(method = 'bfill')
micro_geotype_combined_imputed = pd.concat([micro_geotype_urban_imputed, micro_geotype_rural_imputed])
# micro_geotype_combined_imputed = micro_geotype_combined.fillna(method = 'ffill')
print('total missing after imputation')
print(micro_geotype_combined_imputed.isna().sum())

print('total length after imputation')
print(len(micro_geotype_combined_imputed))

# partition_results_out = partition_results[['GEOID', 'unitedID']]
# partition_results_out['GEOID'] = partition_results_out['GEOID'].astype(str).str.zfill(11)
# micro_geotype_combined_imputed = pd.merge(micro_geotype_combined_imputed,
#                                           partition_results_out, 
#                                           on = 'GEOID', how = 'left')

micro_geotype_combined_imputed.to_csv(os.path.join(path_to_typology, 'microtype_geotype_output_2020.csv'),
                               index = False)


# <codecell>
# impute 2010 assignment
census_tract_2010 = spatial_id_crosswalk_2010[['tract', 'spatial_id']]
census_tract_2010 = census_tract_2010.drop_duplicates(subset = 'tract', keep = 'first')
census_tract_2010['tract'] = \
    census_tract_2010['tract'].astype(str).str.zfill(11)
typology_unimputed_2010['GEOID_TRACT_10'] = \
    typology_unimputed_2010['GEOID_TRACT_10'].astype(str).str.zfill(11)    
micro_geotype_combined_2010 = pd.merge(census_tract_2010, typology_unimputed_2010,
                                       left_on = 'tract', right_on = 'GEOID_TRACT_10',
                                       how = 'left')
print(len(micro_geotype_combined_2010))
micro_geotype_combined_2010 = micro_geotype_combined_2010.rename(columns = {'tract': 'GEOID'})
typology_unimputed_2010 = micro_geotype_combined_2010.sort_values('GEOID', ascending = True)
print('total missing before imputation')
print(typology_unimputed_2010.isna().sum())

#impute rural and urban respectively
micro_geotype_2010_imputed = typology_unimputed_2010.copy()
micro_geotype_2010_imputed.loc[:, 'geotype'] = \
    micro_geotype_2010_imputed.loc[:, 'geotype'].fillna(method = 'ffill')

micro_geotype_urban_imputed = \
    micro_geotype_2010_imputed.loc[micro_geotype_2010_imputed['geotype'].isin(['CBSA_1', 'CBSA_2'])]
micro_geotype_urban_imputed = micro_geotype_urban_imputed.fillna(method = 'ffill')
# micro_geotype_urban_imputed = micro_geotype_urban_imputed.fillna(method = 'bfill')
micro_geotype_rural_imputed = \
    micro_geotype_2010_imputed.loc[~micro_geotype_2010_imputed['geotype'].isin(['CBSA_1', 'CBSA_2'])]
micro_geotype_rural_imputed = micro_geotype_rural_imputed.fillna(method = 'ffill')
# micro_geotype_rural_imputed = micro_geotype_rural_imputed.fillna(method = 'bfill')
micro_geotype_2010_imputed = pd.concat([micro_geotype_urban_imputed, micro_geotype_rural_imputed])
# micro_geotype_combined_imputed = micro_geotype_combined.fillna(method = 'ffill')
print('total missing after imputation')
print(micro_geotype_2010_imputed.isna().sum())

print('total length after imputation')
print(len(micro_geotype_2010_imputed))
# <codecell>      

# APPEND PARTITION TO 2010
#partition_results_out = partition_results_out
census_tract_crosswalk = census_tract_crosswalk[['GEOID_TRACT_20', 'GEOID_TRACT_10']]
census_tract_crosswalk['GEOID_TRACT_20'] = \
    census_tract_crosswalk['GEOID_TRACT_20'].astype(str).str.zfill(11)
census_tract_crosswalk['GEOID_TRACT_10'] = \
    census_tract_crosswalk['GEOID_TRACT_10'].astype(str).str.zfill(11)
    
# partition_results_out_2010 = pd.merge(partition_results_out, census_tract_crosswalk,
#                                        left_on = 'GEOID', right_on = 'GEOID_TRACT_20',
#                                        how = 'left')
# print(len(partition_results_out_2010))

# partition_results_out_2010 = \
#     partition_results_out_2010.drop_duplicates(subset = 'GEOID_TRACT_10', keep = 'first')
# print(len(partition_results_out_2010))   

# # <codecell>      
# partition_results_out_2010 = partition_results_out_2010[['GEOID_TRACT_10', 'unitedID']]         
# micro_geotype_2010_imputed = pd.merge(micro_geotype_2010_imputed,
#                                       partition_results_out_2010,
#                                       on = 'GEOID_TRACT_10',
#                                       how = 'left')     

micro_geotype_2010_imputed = micro_geotype_2010_imputed[['GEOID', 'geotype',
                                                           'network_microtype', 'demand_microtype_comb']]

micro_geotype_2010_imputed.to_csv(os.path.join(path_to_typology, 'microtype_geotype_output_2010.csv'),
                               index = False)         
