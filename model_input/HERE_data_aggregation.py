# -*- coding: utf-8 -*-
"""
Created on Tue Sep  3 11:11:07 2024

@author: xiaodanxu
"""

import pandas as pd
import numpy as np
import os
import matplotlib.pyplot as plt
from pandas import read_csv
import seaborn as sns


data_dir = 'C:\data\CATTLab_delivery'
os.chdir(data_dir)

typology_file = 'AuxiliaryData/microtype_geotype_output_2010.csv'
label_2010 = read_csv(typology_file)

microtype_3_dir = 'National/proc_data/m3_tracts_frc1_5'
microtype_other_dir = 'National/proc_data/all_tracts_frc3_5'

microtype_3_files = [filename for filename in os.listdir(microtype_3_dir) \
                     if filename.endswith("upperbounds.csv")]
microtype_other_files = [filename for filename in os.listdir(microtype_other_dir) \
                     if filename.startswith("bytract_upperbounds_state_")]    
    
label_2010_short = label_2010[['GEOID', 'geotype', 'network_microtype']]
# <codecell>    


HERE_data_microtype_3 = \
    pd.concat((pd.read_csv(os.path.join(microtype_3_dir, f)) for f in microtype_3_files), ignore_index=True)
HERE_data_microtype_3.rename(columns = {'tract_geoid': 'GEOID'}, inplace = True) 
# read microtype 3 data
HERE_data_microtype_3 = pd.merge(HERE_data_microtype_3, label_2010_short, 
                  on = 'GEOID', how = 'left')
# keeping overlapped typology only
HERE_data_microtype_3 = HERE_data_microtype_3.loc[HERE_data_microtype_3['network_microtype'] == 'Urban_2']

def speed_data_processor(df):

    df.loc[:, 'formatted_time'] = \
    pd.to_datetime(df.loc[:, 'datetime'], format="%Y-%m-%d %H:%M:%S")
    df.loc[:, 'weekday'] = df.loc[:, 'formatted_time'].dt.weekday
    df.loc[:, 'hour'] = df.loc[:, 'formatted_time'].dt.hour
    df_out = df.loc[df['weekday'] <= 4]
    df_out.dropna(inplace = True)
    return(df_out)

HERE_data_microtype_3 = speed_data_processor(HERE_data_microtype_3)

 


# <codecell> 
# load data from 5 state
microtype_other_files_sel = ['bytract_upperbounds_state_AL.csv',
                             'bytract_upperbounds_state_CA.csv',
                             'bytract_upperbounds_state_CO.csv',
                             'bytract_upperbounds_state_FL.csv',
                             'bytract_upperbounds_state_ID.csv',
                             'bytract_upperbounds_state_MD.csv']
HERE_data_microtype_other = \
    pd.concat((pd.read_csv(os.path.join(microtype_other_dir, f)) 
               for f in microtype_other_files_sel), ignore_index=True)
HERE_data_microtype_other.rename(columns = {'tract_geoid': 'GEOID'}, inplace = True) 
# read microtype 3 data
HERE_data_microtype_other = pd.merge(HERE_data_microtype_other, 
                                     label_2010_short, 
                  on = 'GEOID', how = 'left')
# keeping overlapped typology only
HERE_data_microtype_other = \
    HERE_data_microtype_other.loc[HERE_data_microtype_other['microtype'] != 3]

HERE_data_microtype_other = \
    HERE_data_microtype_other.loc[HERE_data_microtype_other['network_microtype'] != 'Urban_2']
HERE_data_microtype_other = speed_data_processor(HERE_data_microtype_other)

# <codecell>
HERE_data_by_typology = pd.concat([HERE_data_microtype_3, 
                                   HERE_data_microtype_other])
HERE_data_by_typology.to_csv('Validation/cleaned_HERE_data.csv',
                             index = False)

# <codecell>
plt.style.use('seaborn-v0_8-whitegrid')
# generate hourly avg speed
avg_speed_by_microtype_hour = \
    HERE_data_by_typology.groupby(['geotype', 'network_microtype','hour']).apply(lambda x: np.average(x.aggspeed_mph, weights=x.cnt_tmcs))
avg_speed_by_microtype_hour = avg_speed_by_microtype_hour.reset_index()
avg_speed_by_microtype_hour.columns = ['geotype', 'microtype','hour', 'speed (mph)']
avg_speed_by_microtype_hour['source'] = 'HERE mean speed'
avg_speed_by_microtype_hour.head(5)
avg_speed_by_microtype_hour.to_csv('Validation/avg_speed_from_here_national.csv', index = False)

# <codecell>
# plot hourly mean speed
sns.set_theme(style="whitegrid")
sns.set(font_scale=1.2)  # larger font
# plt.figure(figsize=(8, 6))
sns.relplot(
    data = avg_speed_by_microtype_hour, x = "hour", y = "speed (mph)", 
    col = 'geotype', col_wrap = 2, kind="line", errorbar=None,
    hue="microtype", palette = 'Set1', height = 4, aspect = 1.2,
)
plt.ylim([0, 70])

plt.ylabel('average speed (mph)')

plt.savefig('Validation/national_HERE_hourly_mean_speed.png', 
            bbox_inches='tight', dpi = 300)
plt.show()

# <codecell>

# compute hourly median speed
sns.set_theme(style="white")
plt.style.use('seaborn-v0_8-whitegrid')
def weighted_median(df):
    df_sorted = df.sort_values('aggspeed_mph')
    cumsum = df_sorted['cnt_tmcs'].cumsum()
    cutoff = df_sorted['cnt_tmcs'].sum() / 2.
    median = df_sorted[cumsum >= cutoff]['aggspeed_mph'].iloc[0]
    return(median)

median_speed_by_microtype_hour = \
    HERE_data_by_typology.groupby(['geotype', 'network_microtype','hour']).apply(weighted_median)
median_speed_by_microtype_hour = \
    median_speed_by_microtype_hour.reset_index()
median_speed_by_microtype_hour.columns = ['geotype', 'microtype','hour', 'speed (mph)']
median_speed_by_microtype_hour['source'] = 'HERE median speed'
median_speed_by_microtype_hour.head(5)

median_speed_by_microtype_hour.to_csv('Validation/median_speed_from_here_national.csv', index = False)

# <codecell>
# plot hourly median speed
plt.style.use('seaborn-v0_8-whitegrid')
sns.set_theme(style='white')
sns.set(rc={'axes.facecolor':'white', 'figure.facecolor':'white'})
sns.set(font_scale=1.2)  # larger font
# plt.figure(figsize=(8, 6))
sns.relplot(
    data = median_speed_by_microtype_hour, x = "hour", y = "speed (mph)", 
    col = 'geotype', col_wrap = 2, kind="line", errorbar=None,
    hue="microtype", palette = 'Set1', height = 4, aspect = 1.2,
)
plt.ylim([0, 70])
plt.ylabel('median speed (mph)')
plt.savefig('Validation/national_HERE_hourly_median_speed.png', 
            bbox_inches='tight', dpi = 300)
plt.show()