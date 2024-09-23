# LJ 11/12/2020
# prepare extra variables 
# - percent of trips from a given census under certain distance bin. 
# dist_miles: 
#  <= 2 miles ~ 1, 
#2-5  ~ 2, 
#5-10 ~ 3, 
#10- 20 ~ 4, 
#20 -50 ~ 5, 
#50+ miles ~ 6)

# update 11/17/2020: new distance bins
#mutate(trip_dist_bin = case_when(
#  trpmiles <=1.3  ~ 's', # short
#  trpmiles >1.3 & trpmiles<=3 ~ 'm', # median
#  trpmiles>3 & trpmiles <=8 ~ 'l', # long
#  TRUE ~ 'll')) # very long

#- Trip source magnitude: number of work destinations / number of home origins.
#"this will be the ratio of commute trip destinations over the commute trip origins for each tract,
#- can be derived from: od_pairs_tract_2017.csv.zipped (on the google drive CleanData folder)
#- the variable S100 is the total # of trips we want the ratio of S100 as a destination over the S100 as an origin
##################################################

library(rlist)
library(tidyverse)
library(data.table)

####################################################
# 1. trip distance bins, already computed for each census tract, now need to reshape the data

trips = fread(file.path(inputsdir,'commute_distances_euclidean_tract_updated.csv')) # should be inputsdir

trips = trips %>%
  dplyr::select(tract, dist_bin, pct_trips_bin) %>%
  tidyr::spread(dist_bin,pct_trips_bin)

trips[is.na(trips)]=0 
names(trips) = c('tract',
               'pct_trips_bin1',
               'pct_trips_bin2',
               'pct_trips_bin3',
               'pct_trips_bin4')


# 2. Trip sink magnitude
hw = fread(file.path(inputsdir,'od_pairs_tract_2017.csv')) %>%
  select(w_tract,h_tract,S000)

# now compute for each tract, the ratio of it being a w_tract to it being a h_tract
w = hw %>% 
  dplyr::group_by(w_tract) %>%
  dplyr::summarise(w_tract_total = sum(S000)) %>%
  dplyr::rename(tract = w_tract)

h = hw %>%
  dplyr::group_by(h_tract) %>%
  dplyr::summarise(h_tract_total = sum(S000)) %>%
  dplyr::rename(tract = h_tract)

hw = full_join(w,h) %>%
  filter(!(is.na(h_tract_total)&is.na(w_tract_total))) # remove two tracts that have NA for both h and w

hw$w_tract_total[is.na(hw$w_tract_total)] = 0

hw = hw %>%
  mutate(source_mag = w_tract_total/h_tract_total )

# 3. merge the two sets of additional inputs to stage1 clustering for transgeo paper
add.input = full_join(trips,hw)

#save(add.input,file=file.path(rdatadir,'additional.inputs.RData'))

#Note July 2022: the tract export changed the format in a wierd way, 
#exporting as excel today for the final dataset

write.xlsx(add.input, file.path(inputsdir, "additional.input_2022.xlsx"))


