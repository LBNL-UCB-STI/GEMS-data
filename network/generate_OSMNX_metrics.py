#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Tue Nov 14 10:38:52 2023

@author: xiaodanxu
"""

# from pygris import tracts 
# from pygris import states
import osmnx
import pandas as pd
import geopandas as gpd
from multiprocessing import Pool
import time
import os

def split_dataframe(df, chunk_size = 10 ** 6): 
    chunks = list()
    num_chunks = len(df) // chunk_size + 1
    for i in range(num_chunks):
        chunks.append(df[i*chunk_size:(i+1)*chunk_size])
    return chunks

def query_osmnx_stats(geometry, var_to_keep):
    try:
        G = osmnx.graph_from_polygon(geometry, network_type='all')
        # osmnx.plot_graph(G)
        basic_stats = pd.Series(osmnx.basic_stats(G))
        basic_stats = basic_stats[var_to_keep]
    except:
        basic_stats = pd.Series(index=var_to_keep)
    return(basic_stats)

def add_df_attr(args):
    df, i, var_to_keep = args
    out_path = 'Network/RawData/OSMNX/osmnx_stats_by_tract_chunk_' + str(i) + '.csv'
    if os.path.exists(out_path):
        return True
    df.loc[:, var_to_keep] = \
        df.apply(lambda row : query_osmnx_stats(row['geometry'], var_to_keep), axis = 1)
    df = df.drop(columns = 'geometry')
    df.to_csv(out_path, index = False)
    return(df)

def main():
    print('run start')
    path = 'C:/FHWA_R2'
    analysis_year = 2020
    os.chdir(path)
    start_time = time.time()
    # us_states = states(cb=True, year=2022)
    # state_to_drop = ['AS', 'PR', 'GU', 'VI', 'MP']
    # us_states = us_states.loc[~us_states['STUSPS'].isin(state_to_drop)]
    # list_of_state = us_states['STUSPS'].unique()
    # state_tracts = None
    # for st in list_of_state:
    #     state_tracts_sel = tracts(state = st, cb = False, year = 2022, cache = True)
    #     state_tracts = pd.concat([state_tracts, state_tracts_sel])
    # print(len(state_tracts))
    state_tracts = gpd.read_file('spatial_boundary/CleanData/combined_tracts_' + str(analysis_year) +'.geojson')
    print('finish map loading')
    # <codecell>
    
    chunk_size = 10 ** 3
    njob = 0
    var_to_keep = ['n', 'm', 'k_avg', 'edge_length_total', 'edge_length_avg', 
                   'streets_per_node_avg', 'intersection_count', 'street_length_total', 
                   'street_segment_count', 
                   'street_length_avg', 'circuity_avg', 'self_loop_proportion']
    # sample_tracts_test = state_tracts.head(10)
    chunks_of_tracts = split_dataframe(state_tracts, chunk_size)
    # CA_tracts_test.plot()
    i = 0
    jobs=[]
    for chunk in chunks_of_tracts:
        jobs.append( (chunk, i, var_to_keep))
        i+=1

    print(len(jobs))
    print(jobs[0])
    njob+=len(jobs)
    pl=Pool(8)
    pl.map(add_df_attr, jobs)
    
    # sample_tracts_test.loc[:, var_to_keep] = \
    #     sample_tracts_test.apply(lambda row : query_osmnx_stats(row['geometry'], var_to_keep), axis = 1)
    end_time = time.time()
    
    print('total computation time')
    print(end_time - start_time)
# G_feature = osmnx.features_from_polygon(CA_tracts_test.iloc[0]['geometry'], 
#                                         tags = {'highway': True})

if __name__ == '__main__':
	main()