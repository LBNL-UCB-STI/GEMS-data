# -*- coding: utf-8 -*-
"""
Created on Thu Nov 16 10:55:38 2023

@author: xiaodanxu
"""

from pandas import read_csv
import os
from os import listdir
import pandas as pd

var_to_keep = ['st', 'stusps', 'stname', 'cty', 'ctyname', 'trct', 'trctname', 'cbsa', 'cbsaname']
os.chdir('C:/FHWA_R2/spatial_boundary')
all_files = os.listdir("RawData/LODES8") 
csv_files = list(filter(lambda f: f.endswith('.csv'), all_files))
cleaned_crosswalk = None
for file in csv_files:
    print('processing file ' + file)
    df = read_csv("RawData/LODES8/" + file) 
    df = df[var_to_keep]
    # print(len(df))
    df = df.drop_duplicates(keep = 'first')
    # print(len(df))
    cleaned_crosswalk = pd.concat([cleaned_crosswalk, df])

# <codecell>
cleaned_crosswalk.to_csv('CleanData/cleaned_lodes8_crosswalk.csv', index = False)
