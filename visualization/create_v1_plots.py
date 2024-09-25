#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Mon Sep 16 13:07:52 2024

@author: xiaodanxu
"""

import pandas as pd
from pandas import read_csv
import numpy as np
import matplotlib.pyplot as plt
import seaborn as sns
import os
import shapely.wkt
import geopandas as gpd
import contextily as cx
import warnings
warnings.filterwarnings("ignore")


# plot V1 typology
v1_cluster = read_csv('data/ccst_geoid_key_transp_geo_with_imputation.csv')
v1_cluster_result_count = v1_cluster.groupby(['geotype', 'microtype']).size()
# v1_cluster_result_count.head(5)

v1_cluster_result_count = v1_cluster_result_count.reset_index()
v1_cluster_result_count.columns = ['geotype', 'microtype', 'count']
v1_cluster_result_count.head(5)

sns.set_theme(style="white", font_scale=1.4)
ax = sns.catplot(data = v1_cluster_result_count, 
            x="microtype", y="count", 
            col="geotype", 
            col_wrap = 3, kind="bar", 
            # hue_order = order_network,
            # order = order_demand,
            palette = 'YlOrRd_r',
            height = 5, aspect = 1.4, sharey = True)
# ax.set_titles("{col_name}")
# plt.xticks(rotation = 90, ha= 'right')
for axn in ax.axes.flat:
    for label in axn.get_xticklabels():
        label.set_rotation(0)
plt.savefig('plot/V1_tract_by_typology.png', dpi = 300,
           bbox_inches = 'tight')
plt.show()
