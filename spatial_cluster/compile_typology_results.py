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

network_typology_cbsa_file = 'network_geotype_cbsa_clustered_exEdgeCount.csv'
network_typology_noncbsa_file = 'network_geotype_non_cbsa_clustered_exEdgeCount.csv'
partition_results_file = 'final_partition_results.csv'
urban_geotype_file = 'Geotype_Urban_Clusters_k2_4_pam.csv'
rural_geotype_file = 'Geotype_Rural_Clusters_k2_4_pam.csv'
demand_microtype_file = 'clustering_outputs_with_raw_data.csv'

path_to_crosswalk = 'spatial_boundary/CleanData/'
census_tract_crosswalk_file = 'census_tract_crosswalk_2010_2020.csv'
spatial_id_crosswalk_file = 'cleaned_lodes8_crosswalk_with_ID.csv'

network_typology_cbsa = read_csv(os.path.join(path_to_typology, network_typology_cbsa_file), low_memory=False)
network_typology_noncbsa = read_csv(os.path.join(path_to_typology, network_typology_noncbsa_file), low_memory=False)
partition_results = read_csv(os.path.join(path_to_typology, partition_results_file))

urban_geotype = read_csv(os.path.join(path_to_typology, urban_geotype_file))
rural_geotype = read_csv(os.path.join(path_to_typology, rural_geotype_file))

demand_microtype = read_csv(os.path.join(path_to_demand, demand_microtype_file))

census_tract_crosswalk = read_csv(os.path.join(path_to_crosswalk, census_tract_crosswalk_file))
spatial_id_crosswalk = read_csv(os.path.join(path_to_crosswalk, spatial_id_crosswalk_file))

network_urban_clusters = 'cluster5'
network_rural_clusters = 'cluster2'

geotype_urban_clusters = 'cluster2'
geotype_rural_clusters = 'cluster2'
# <codecell>

# matching urban network typology to tract
partition_ids = partition_results.unitedID.unique()

network_typology_cbsa = network_typology_cbsa[['network_id', network_urban_clusters]]

network_typology_cbsa_par = \
    network_typology_cbsa.loc[network_typology_cbsa['network_id'].isin(partition_ids)]
    
network_typology_cbsa_nopar = \
    network_typology_cbsa.loc[~network_typology_cbsa['network_id'].isin(partition_ids)]  

network_typology_cbsa_par = pd.merge(network_typology_cbsa_par,
                                     partition_results, left_on = 'network_id', 
                                     right_on = 'unitedID', how = 'left')

network_typology_cbsa_par = network_typology_cbsa_par[['GEOID', network_urban_clusters]]

network_typology_cbsa_nopar = network_typology_cbsa_nopar.rename(columns = {'network_id': 'GEOID'})

network_typology_cbsa_processed = pd.concat([network_typology_cbsa_par,
                                             network_typology_cbsa_nopar])
# <codecell>
# appending rural network typology

print(network_typology_noncbsa.columns)

network_typology_noncbsa = network_typology_noncbsa[['network_id', network_rural_clusters]]

network_typology_noncbsa = network_typology_noncbsa.rename(columns = {'network_id': 'GEOID'})

network_typology_cbsa_processed.loc[:, 'network_microtype'] = 'Urban_' + \
    network_typology_cbsa_processed.loc[:, network_urban_clusters].astype(str)
    
network_typology_noncbsa.loc[:, 'network_microtype'] = 'Rural_' + \
    network_typology_noncbsa.loc[:, network_rural_clusters].astype(str)

output_attr = ['GEOID', 'network_microtype']
network_typology_output = pd.concat([network_typology_cbsa_processed[output_attr],
                                     network_typology_noncbsa[output_attr]])
print(len(network_typology_output)) # 83612
network_typology_output.to_csv(os.path.join(path_to_typology, 'network_typology_output_2020.csv'),
                               index = False)

# <codecell>

# assign geotype to census tract

print(urban_geotype.columns)
urban_geotype = urban_geotype[['spatial_id', geotype_urban_clusters]]
rural_geotype = rural_geotype[['spatial_id', geotype_rural_clusters]]

urban_geotype.loc[:, 'geotype'] = 'CBSA_' + \
    urban_geotype.loc[:, geotype_urban_clusters].astype(str)
    
rural_geotype.loc[:, 'geotype'] = 'NONCBSA_' + \
    rural_geotype.loc[:, geotype_rural_clusters].astype(str)

output_attr = ['spatial_id', 'geotype']
geotype_combined = pd.concat([urban_geotype[output_attr], 
                              rural_geotype[output_attr]])

spatial_id_crosswalk = spatial_id_crosswalk[['spatial_id', 'trct']]

geotype_combined = pd.merge(spatial_id_crosswalk, geotype_combined, 
                            on = 'spatial_id', how = 'left')

geotype_combined = geotype_combined.rename(columns = {'trct':'GEOID'})
print(len(geotype_combined)) # 84412
geotype_combined.to_csv(os.path.join(path_to_typology, 'geotype_output_2020.csv'),
                               index = False)

# <codecell>

# consolidate micro-geotype and convert to 2010 boundary
demand_microtype = demand_microtype[['GEOID', 'demand_microtype_comb']]
geotype_combined['GEOID'] = geotype_combined['GEOID'].astype(str).str.zfill(11)
network_typology_output['GEOID'] = network_typology_output['GEOID'].astype(str).str.zfill(11)
demand_microtype['GEOID'] = demand_microtype['GEOID'].astype(str).str.zfill(11)

micro_geotype_combined = pd.merge(geotype_combined, demand_microtype, 
                                  on = 'GEOID', how = 'outer')

micro_geotype_combined = pd.merge(micro_geotype_combined, network_typology_output,
                                  on = 'GEOID', how = 'outer')

micro_geotype_combined.loc[micro_geotype_combined['spatial_id'] == 15005, 'geotype'] = 'County_1'
print(len(micro_geotype_combined))
micro_geotype_combined.to_csv(os.path.join(path_to_typology, 'microtype_geotype_output_2020_unimputed.csv'),
                               index = False)

# <codecell>
# fill missing value
micro_geotype_combined = micro_geotype_combined.sort_values('GEOID', ascending = True)
print('total missing before imputation')
print(micro_geotype_combined.isna().sum())

#impute rural and urban respectively
micro_geotype_combined_imputed = micro_geotype_combined.copy()
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

micro_geotype_combined_imputed.to_csv(os.path.join(path_to_typology, 'microtype_geotype_output_2020.csv'),
                               index = False)
# <codecell>
census_tract_crosswalk = census_tract_crosswalk[['GEOID_TRACT_20', 'GEOID_TRACT_10']]
census_tract_crosswalk['GEOID_TRACT_20'] = \
    census_tract_crosswalk['GEOID_TRACT_20'].astype(str).str.zfill(11)
census_tract_crosswalk['GEOID_TRACT_10'] = \
    census_tract_crosswalk['GEOID_TRACT_10'].astype(str).str.zfill(11)
    
micro_geotype_combined_2010 = pd.merge(micro_geotype_combined_imputed, census_tract_crosswalk,
                                       left_on = 'GEOID', right_on = 'GEOID_TRACT_20',
                                       how = 'left')
print(len(micro_geotype_combined_2010))
# <codecell>
micro_geotype_combined_2010 = \
    micro_geotype_combined_2010.drop_duplicates(subset = 'GEOID_TRACT_10', keep = 'first')
print(len(micro_geotype_combined_2010))

micro_geotype_combined_2010 = micro_geotype_combined_2010[['GEOID_TRACT_10', 'geotype',
                                                           'network_microtype', 'demand_microtype_comb']]

micro_geotype_combined_2010.to_csv(os.path.join(path_to_typology, 'microtype_geotype_output_2010.csv'),
                               index = False)






