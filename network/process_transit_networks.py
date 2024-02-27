# -*- coding: utf-8 -*-
"""
Created on Fri Feb 16 10:15:17 2024
for processing NTAD transit data to census tract level
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
# NTAD transit data accessed Feb 16, 2024
stop_file = os.path.join('Network', data_path, 'NTAD', 'National_Transit_Map_Stops.geojson')
national_transit_stop = gpd.read_file(stop_file)

route_file = os.path.join('Network', data_path, 'NTAD', 'National_Transit_Map_Routes.geojson')
national_transit_route = gpd.read_file(route_file)

# <codecell>

# load census tract boundary

# use 2018 boundary for mode choice estimation and 2020 for GEMS 2.0
run_type_desc = {1: 'mode choice', 2: 'GEMS update'}
run_type = 1
if run_type == 1:
    ct_file = 'spatial_boundary/CleanData/combined_tracts_2018.geojson'
else:
    ct_file = 'spatial_boundary/CleanData/combined_tracts_2020.geojson'
    
us_census_tract = gpd.read_file(ct_file)
# change projection
ntad_crs = national_transit_route.crs
us_census_tract = us_census_tract.to_crs(ntad_crs)

# <codecell>

#bus route selection
bus_types = ['3', '11', '5']
national_bus_route = national_transit_route.loc[national_transit_route['route_type'].isin(bus_types)]

rail_types = ['0', '1', '2', '12']
national_rail_route = national_transit_route.loc[national_transit_route['route_type'].isin(rail_types)]

# <codecell>

# assign rail to census tract


rail_to_tract = gpd.overlay(national_rail_route, us_census_tract, how='intersection')

# compute segment length in meters
rail_to_tract = rail_to_tract.to_crs("EPSG:3310") 
# in order to get length in meter, the shapefile need to re-projected to a coordinate system in meters (not required in R)
rail_to_tract.loc[:, 'Length'] = rail_to_tract.loc[:, 'geometry'].length

# <codecell>
rail_to_tract_df = pd.DataFrame(rail_to_tract.drop(columns='geometry'))

# create census tract DF
us_census_tract_df =  pd.DataFrame(us_census_tract.drop(columns='geometry'))
rail_dist_by_tract = rail_to_tract_df.groupby(['GEOID'])[['Length']].sum()
rail_dist_by_tract.columns = ['rail_length_m']
rail_dist_by_tract = rail_dist_by_tract.reset_index()
us_census_tract_df = pd.merge(us_census_tract_df, rail_dist_by_tract,
                              on = 'GEOID', how = 'left')

# <codecell>
ca_tract = pygris.tracts(state = "CA", year = 2018)

ca_tract = ca_tract.merge(rail_dist_by_tract, 
                          on = 'GEOID', how ='left')
ca_tract.to_file(os.path.join('Network', output_path, "sample_ca_rail_dist.geojson"), driver='GeoJSON')

rail_to_tract_ca = rail_to_tract.loc[rail_to_tract['STATEFP'] == '06']
rail_to_tract_ca.to_file(os.path.join('Network', output_path, "sample_ca_rail_route.geojson"), driver='GeoJSON')

# <codecell>
# assign bus to census tract

bus_to_tract = gpd.overlay(national_bus_route, us_census_tract, how='intersection')

# compute segment length in meters
bus_to_tract = bus_to_tract.to_crs("EPSG:3310") 
# in order to get length in meter, the shapefile need to re-projected to a coordinate system in meters (not required in R)
bus_to_tract.loc[:, 'Length'] = bus_to_tract.loc[:, 'geometry'].length

# <codecell>
bus_to_tract_df = pd.DataFrame(bus_to_tract.drop(columns='geometry'))

# create census tract DF
# us_census_tract_df =  pd.DataFrame(us_census_tract.drop(columns='geometry'))
bus_dist_by_tract = bus_to_tract_df.groupby(['GEOID'])[['Length']].sum()
bus_dist_by_tract.columns = ['bus_length_m']
bus_dist_by_tract = bus_dist_by_tract.reset_index()
us_census_tract_df = pd.merge(us_census_tract_df, bus_dist_by_tract,
                              on = 'GEOID', how = 'left')

# <codecell>
ca_tract = pygris.tracts(state = "CA", year = 2018)

ca_tract = ca_tract.merge(bus_dist_by_tract, 
                          on = 'GEOID', how ='left')
ca_tract.to_file(os.path.join('Network', output_path, "sample_ca_bus_dist.geojson"), driver='GeoJSON')

bus_to_tract_ca = bus_to_tract.loc[bus_to_tract['STATEFP'] == '06']
bus_to_tract_ca.to_file(os.path.join('Network', output_path, "sample_ca_bus_route.geojson"), driver='GeoJSON')

# <codecell>

# clean data and write output
us_census_tract_df = us_census_tract_df.fillna(0)
us_census_tract_df.loc[:, 'rail_density'] = \
    us_census_tract_df.loc[:, 'rail_length_m'] / us_census_tract_df.loc[:, 'ALAND']
    
us_census_tract_df.loc[:, 'bus_density'] = \
    us_census_tract_df.loc[:, 'bus_length_m'] / us_census_tract_df.loc[:, 'ALAND']
    
us_census_tract_df.loc[:, 'rail_available'] = 0
us_census_tract_df.loc[us_census_tract_df['rail_length_m'] > 0, 'rail_available'] = 1

us_census_tract_df.loc[:, 'bus_available'] = 0
us_census_tract_df.loc[us_census_tract_df['bus_length_m'] > 0, 'bus_available'] = 1

us_census_tract_df.to_csv(os.path.join('Network', output_path, "transit_availability.csv"), index = False)


# <codecell>

# compare with V1 
v1_file = os.path.join('Network', data_path, 'NTAD', 'modeaccessbility.csv')
mode_availability_v1 = read_csv(v1_file)
us_census_tract_df['GEOID'] = us_census_tract_df['GEOID'].astype(str).str.zfill(11)
mode_availability_v1['geoid'] = mode_availability_v1['geoid'].astype(str).str.zfill(11)

us_census_tract_compare = pd.merge(us_census_tract_df, mode_availability_v1,
                                   left_on = 'GEOID', right_on = 'geoid', how = 'left')

print(us_census_tract_compare.head(5))

from sklearn.metrics import confusion_matrix
rail_cross_tab = confusion_matrix( us_census_tract_compare['rail'], 
                                  us_census_tract_compare['rail_available'])
print(rail_cross_tab)
bus_cross_tab = confusion_matrix( us_census_tract_compare['bus'],
                                 us_census_tract_compare['bus_available'])
print(bus_cross_tab)

# <codecell>

# <changed the way to calculate availability
us_census_tract = us_census_tract.to_crs("EPSG:3310") 
us_census_tract_centroid = us_census_tract.centroid
us_census_tract_centroid = pd.concat([us_census_tract_df, us_census_tract_centroid], axis = 1)
us_census_tract_centroid = us_census_tract_centroid.rename(columns = {0: 'geometry'})
us_census_tract_centroid = gpd.GeoDataFrame(us_census_tract_centroid, geometry='geometry')

# <codecell>
# calculate nearest dist to rail
national_rail_route = national_rail_route.to_crs("EPSG:3310")
us_census_tract_centroid.loc[:, 'dist_to_rail_m'] = \
    us_census_tract_centroid.geometry.apply(lambda x: national_rail_route.distance(x).min())

# <codecell>
us_census_tract_centroid.loc[:, 'dist_to_rail_m'] = us_census_tract_centroid.loc[:, 'dist_to_rail_m'].astype(float)
# us_census_tract_centroid[:, 'dist_to_rail_mile'] = 0.000621371 * us_census_tract_centroid.loc[:, 'dist_to_rail_m']
us_census_tract_centroid_df = pd.DataFrame(us_census_tract_centroid.drop(columns='geometry'))
us_census_tract_centroid_df.to_csv(os.path.join('Network', output_path, "transit_availability_with_dist.csv"), index = False)

    # calculate nearest dist to bus
# national_bus_route = national_bus_route.to_crs("EPSG:3310")
# us_census_tract_centroid.loc[:, 'dist_to_bus_m'] = \
#     us_census_tract_centroid.geometry.apply(lambda x: national_bus_route.distance(x).min())


