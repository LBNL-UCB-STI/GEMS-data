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

data_aggregation_done = 0 # 1- yes, 0 - no

microtype_3_files = [filename for filename in os.listdir(microtype_3_dir) \
                     if filename.endswith("upperbounds.csv")]
microtype_other_files = [filename for filename in os.listdir(microtype_other_dir) \
                     if filename.startswith("bytract_upperbounds_state_")]    
    
label_2010_short = label_2010[['GEOID', 'geotype', 'network_microtype']]
# <codecell>    
meter_per_mile = 1609.34
def speed_data_processor(df):

    df.loc[:, 'formatted_time'] = \
    pd.to_datetime(df.loc[:, 'datetime'], format="%Y-%m-%d %H:%M:%S")
    df.loc[:, 'weekday'] = df.loc[:, 'formatted_time'].dt.weekday
    df.loc[:, 'hour'] = df.loc[:, 'formatted_time'].dt.hour
    df_out = df.loc[df['weekday'] <= 4]
    thres_density = 100 / meter_per_mile
    df_out = df_out.loc[df_out['agg_density_plpm'] >= thres_density]
    df_out.dropna(inplace = True)
    return(df_out)

def weighted_median(df):
    df_sorted = df.sort_values('aggspeed_mph')
    cumsum = df_sorted['cnt_tmcs'].cumsum()
    cutoff = df_sorted['cnt_tmcs'].sum() / 2.
    median = df_sorted[cumsum >= cutoff]['aggspeed_mph'].iloc[0]
    return(median)

if data_aggregation_done == 0:
    HERE_data_microtype_3 = \
        pd.concat((pd.read_csv(os.path.join(microtype_3_dir, f)) for f in microtype_3_files), ignore_index=True)
    HERE_data_microtype_3.rename(columns = {'tract_geoid': 'GEOID'}, inplace = True) 
    # read microtype 3 data
    HERE_data_microtype_3 = pd.merge(HERE_data_microtype_3, label_2010_short, 
                      on = 'GEOID', how = 'left')
    # keeping overlapped typology only
    HERE_data_microtype_3 = HERE_data_microtype_3.loc[HERE_data_microtype_3['network_microtype'] == 'Urban_2']
    
    
    
    HERE_data_microtype_3 = speed_data_processor(HERE_data_microtype_3)

    # load data from 5 state
    microtype_other_files_sel = ['bytract_upperbounds_state_AL.csv',
                                 'bytract_upperbounds_state_CA.csv',
                                 'bytract_upperbounds_state_CO.csv',
                                 'bytract_upperbounds_state_CT.csv',
                                 'bytract_upperbounds_state_FL.csv',
                                 'bytract_upperbounds_state_ID.csv',
                                 'bytract_upperbounds_state_MD.csv',
                                 'bytract_upperbounds_state_KY.csv',
                                 'bytract_upperbounds_state_SD.csv']
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
    
    HERE_data_by_typology = pd.concat([HERE_data_microtype_3, 
                                       HERE_data_microtype_other])
    HERE_data_by_typology.to_csv('Validation/cleaned_HERE_data.csv',
                                 index = False)
    # generate hourly avg speed
    avg_speed_by_microtype_hour = \
        HERE_data_by_typology.groupby(['geotype', 'network_microtype','hour']).apply(lambda x: np.average(x.aggspeed_mph, weights=x.cnt_tmcs))
    avg_speed_by_microtype_hour = avg_speed_by_microtype_hour.reset_index()
    avg_speed_by_microtype_hour.columns = ['geotype', 'microtype','hour', 'speed (mph)']
    avg_speed_by_microtype_hour['source'] = 'HERE mean speed'
    avg_speed_by_microtype_hour.head(5)
    avg_speed_by_microtype_hour.to_csv('Validation/avg_speed_from_here_national.csv', index = False)
    # compute hourly median speed
    median_speed_by_microtype_hour = \
        HERE_data_by_typology.groupby(['geotype', 'network_microtype','hour']).apply(weighted_median)
    median_speed_by_microtype_hour = \
        median_speed_by_microtype_hour.reset_index()
    median_speed_by_microtype_hour.columns = ['geotype', 'microtype','hour', 'speed (mph)']
    median_speed_by_microtype_hour['source'] = 'HERE median speed'
    median_speed_by_microtype_hour.head(5)

    median_speed_by_microtype_hour.to_csv('Validation/median_speed_from_here_national.csv', index = False)
    
else:
    avg_speed_by_microtype_hour = read_csv('Validation/avg_speed_from_here_national.csv')
    median_speed_by_microtype_hour = read_csv('Validation/median_speed_from_here_national.csv')
# <codecell>
plt.style.use('seaborn-v0_8-whitegrid')


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

# <codecell>

# load GEMS data
mps_to_mph = 2.23694
list_of_gems_output_dir = ['Validation/newdata-all-geotypes-CBSA_1/', 
                           'Validation/newdata-all-geotypes-CBSA_2/',
                          'Validation/newdata-all-geotypes-NONCBSA_1/',
                          'Validation/newdata-all-geotypes-NONCBSA_2/']
speed_results_file = 'speed.csv'

combined_GEMS_speed_output = None
for gems_output_dir in list_of_gems_output_dir:
    file_dir = gems_output_dir.split('/')[1]
    geotype = file_dir.split('-')[3]
    print(geotype)
    GEMS_speed_output = pd.read_csv(gems_output_dir + speed_results_file, header = [0])
    # GEMS_speed_output.columns = GEMS_speed_output.columns.map('_'.join).str.strip('_')
    current_names = GEMS_speed_output.columns
    GEMS_speed_output = GEMS_speed_output.rename(columns = {current_names[0]: 'microtype', current_names[1]: 'mode'})
    GEMS_speed_output.loc[:, 'geotype'] = geotype
    combined_GEMS_speed_output = pd.concat([combined_GEMS_speed_output, GEMS_speed_output])
print(combined_GEMS_speed_output.columns)

# <codecell>

# load time bin data
param_dir = 'AuxiliaryData/'
time_period_file = 'TimePeriods.csv'
nSubBins = 1

time_period_spec = pd.read_csv(param_dir + time_period_file, sep = ',')
time_period_spec.loc[:, 'Interval'] = time_period_spec.loc[:, 'DurationInHours'] / nSubBins
time_period_spec.loc[:, 'EndTime']  = time_period_spec.loc[:, 'DurationInHours'].cumsum()
gen_cols = []
for i in range(nSubBins):
    #print(i)
    colname = 't' + str(nSubBins - i)
    gen_cols.append(colname)
    time_period_spec.loc[:, colname] = time_period_spec.loc[:, 'EndTime'] - (i+1) * time_period_spec.loc[:, 'Interval'] 

time_period_spec = pd.melt(time_period_spec, id_vars=['TimePeriodID'], value_vars=gen_cols, value_name = 'hour')
time_period_spec = time_period_spec.sort_values('hour')
time_period_spec = time_period_spec.loc[:, ['hour']]
time_period_spec.loc[:, 'ID'] = np.arange(time_period_spec.shape[0])
print(time_period_spec.head(5))

# <codecell>

renaming_network_type = {'Urban_1': 'Urban_5',
                         'Urban_2': 'Urban_4',
                         'Urban_3': 'Urban_3',
                         'Urban_4': 'Urban_2',
                         'Urban_5': 'Urban_1',
                         'Rural_1': 'Rural_3',
                         'Rural_2': 'Rural_2',
                         'Rural_3': 'Rural_1'
    }

renaming_geotype = {'CBSA_1': 'B',
                    'CBSA_2': 'A',
                    'NONCBSA_1': 'C',
                    'NONCBSA_2': 'D'}

combined_GEMS_speed_output.loc[:, 'geotype'] = \
    combined_GEMS_speed_output.loc[:, 'microtype'].str.split('-').str[0]
    
combined_GEMS_speed_output.loc[:, 'geotype'] = \
    combined_GEMS_speed_output.loc[:, 'geotype'].map(renaming_geotype)
    
combined_GEMS_speed_output.loc[:, 'microtype'] = \
    combined_GEMS_speed_output.loc[:, 'microtype'].str.split('-').str[1]
    
combined_GEMS_speed_output.loc[:, 'microtype'] = \
    combined_GEMS_speed_output.loc[:, 'microtype'].map(renaming_network_type)
    
print(combined_GEMS_speed_output.head(5))

# <codecell>
list_of_speed = ['Speed_' + str(i)  for i in np.arange(time_period_spec.shape[0])]
# print(list_of_speed)
mode_to_select = 'auto'
GEMS_speed_output_auto = combined_GEMS_speed_output.loc[combined_GEMS_speed_output['mode'] == mode_to_select]
GEMS_speed_output_auto = pd.melt(GEMS_speed_output_auto, id_vars = ['geotype', 'microtype'], 
                                 value_vars = list_of_speed, value_name = 'speed')
GEMS_speed_output_auto = GEMS_speed_output_auto.reset_index()
GEMS_speed_output_auto.loc[:, 'speed'] *= mps_to_mph
GEMS_speed_output_auto.loc[:, 'ID'] = GEMS_speed_output_auto.loc[:, 'variable'].str.split('_').str[1]
GEMS_speed_output_auto.loc[:, 'ID'] = GEMS_speed_output_auto.loc[:, 'ID'].astype(int)
GEMS_speed_output_auto = pd.merge(GEMS_speed_output_auto, time_period_spec, on = ['ID'], how = 'left')
# print(GEMS_speed_output_auto.head(5))

GEMS_speed_output_auto_to_plot = GEMS_speed_output_auto.loc[:, ['geotype', 'microtype', 'hour', 'speed']]
GEMS_speed_output_auto_to_plot.columns = ['geotype', 'network type', 'hour', 'speed (mph)']
GEMS_speed_output_auto_to_plot.loc[:,'source'] = 'GEMS lower-level'

median_speed_by_microtype_hour.rename(columns = {'microtype': 'network type'}, inplace = True)
avg_speed_by_microtype_hour.rename(columns = {'microtype': 'network type'}, inplace = True)
speed_to_plot = pd.concat([median_speed_by_microtype_hour, avg_speed_by_microtype_hour, GEMS_speed_output_auto_to_plot])
speed_to_plot = speed_to_plot.reset_index()
print(speed_to_plot.head(5))

# <codecell>
list_of_geotype = GEMS_speed_output_auto_to_plot.geotype.unique()
for gt in list_of_geotype:
    speed_to_plot_by_gt = speed_to_plot.loc[speed_to_plot['geotype'] == gt]
    ax = sns.relplot(x="hour", y="speed (mph)", hue="source", style="source", col="network type", col_wrap = 2,
        height=2.5, aspect=2, kind="line", linewidth=2.5, palette = 'tab20', data=speed_to_plot_by_gt)
    plt.ylim([0, 70])
    # plt.ylabel('Speed (mph)')
    # plt.title('Sample speed comparison')
    # plt.legend(loc='center right', title='Source')
    # ax.legend(loc='center left', bbox_to_anchor=(1, 0.5))
    plt.savefig('Validation/plot/' + 'speed_comparison_national_geotype_' + gt + '.png', bbox_inches='tight', dpi = 300)
    plt.show()
    
# <codecell>
def time_series_smape(df):
    observed_speed = df.loc[df['source'] == 'HERE mean speed', 'speed (mph)'].array
    modeled_speed = df.loc[df['source'] == 'GEMS lower-level', 'speed (mph)'].array
    ratio = abs(observed_speed - modeled_speed)/(abs(modeled_speed) + abs(observed_speed))
    MAPE = np.mean(ratio)
    return(MAPE)

# speed_to_plot_geotype_A = speed_to_plot.loc[speed_to_plot['geotype'] == 'A']
for gt in list_of_geotype:
    speed_to_plot_selected = speed_to_plot.loc[speed_to_plot['geotype'] == gt]
    time_series_mape = speed_to_plot_selected.groupby(['network type']).apply(time_series_smape)
    time_series_mape = time_series_mape.reset_index()
    time_series_mape.columns = ['network type', 'time series MAPE']
    time_series_mape.to_csv('Validation/plot/' + 'national_smape_by_geotype' + gt + '.csv', index = False)

time_series_mape.head(6)
