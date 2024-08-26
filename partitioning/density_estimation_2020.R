# Rscripts to test the spatial coverage by microtypes in 7 cities.

# set working directory
mywd <- "D:/Projects/Berkeley/GEMS/" 
setwd(mywd)
library(pacman)
p_load(R.utils) # for reading gz files
p_load(readr) # for read files from a zip foler

p_load(data.table)
p_load(dplyr)
p_load(tidyr)
p_load(ggplot2)
p_load(stringr)
p_load(vroom)
p_load(tictoc)
p_load(parallel)
p_load(doParallel)
p_load(foreach)
p_load(stargazer)
p_load(Hmisc) # weighted quantile
p_load(stringr)
p_load(readr) # read multiple csv
p_load(bit64) # convert char to integer64
p_load(sf)
p_load(rgeoda) # for spatial clustering
p_load(geodaData) # sample data for spatial clustering
p_load(units) # to drop units
# p_load(rmapshaper) # to remove isolated polygons, not working. so I write my own code to do it
# #p_load(geojson_sf) # convert to geojson and sf, but this is not available
p_load(tmap) # visualize spatial features
p_load(clusterCrit) # validate cluster
p_load(tigris)
p_load(rgdal)
p_load(raster)
p_load(ggmap)
p_load(maptiles)
p_load(tidyterra)
p_load(hexbin) # viridis color
p_load(viridisLite)
p_load(grid)
p_load(mlogit)
p_load(gridExtra) # plot multiple plots
source('functions.R')
source('spatial.partition_functions.R')

# set directories
auxdatdir <- 'D:/Projects/Berkeley/GEMS/auxiliary_data/' # the blob-microtype ids
density20 <- 'G:/My Drive/MFD/Tract_level_clustering/National/proc_data/density'
rdatdir <-'G:/My Drive/MFD/Tract_level_clustering/National/proc_data' # where you put your final results

# if on Ling's Mac: this is where to find the flow density data
# otherwise, set it to your directory where you downloaded the data
cat.proc <- 'G:/My Drive/MFD/Tract_level_clustering/National/proc_data/tract_mr16'

df <- read_csv("D:/Projects/Berkeley/GEMS/typology/network_cluster_input_updated.csv")
df2 <- read_csv("D:/Projects/Berkeley/GEMS/typology/census_tract_crosswalk_2010_2020.csv")
df3 <- df2[ ,c('GEOID_TRACT_10', 'AREALAND_TRACT_10')]
df3 <- unique(df3)

extracted_geoids <- df[df$cbsa != 99999, 'GEOID']
subset_df2 <- df2[df2$GEOID_TRACT_20 %in% extracted_geoids$GEOID, c('GEOID_TRACT_20', 'GEOID_TRACT_10', 'AREALAND_TRACT_10')]

# load microtype definition by tract
linker = fread(file.path(auxdatdir,'ccst_geoid_key_tranps_geo_with_imputation.csv')) %>%
  dplyr::select(GEOID, st_code)%>%
  distinct()

# get weights for matching GEOID2010
weights_df <- subset_df2 %>%
  group_by(GEOID_TRACT_20) %>%
  mutate(Total_AREALAND = sum(AREALAND_TRACT_10)) %>%
  ungroup() %>%
  group_by(GEOID_TRACT_20, GEOID_TRACT_10) %>%
  summarize(Weight = first(AREALAND_TRACT_10) / first(Total_AREALAND), .groups = 'drop') %>%
  group_by(GEOID_TRACT_20) %>%
  summarise(GEOID_TRACT_10 = list(GEOID_TRACT_10),
            Weights = list(Weight))


contains_nan <- function(x) {
  if (is.numeric(x)) {
    return(any(is.nan(x)))  # Returns TRUE if any value is NaN
  } else {
    return(FALSE)  # Returns FALSE for non-numeric data (adjust based on your data characteristics)
  }
}

# Apply the function to each element of the list column to identify rows with NaN
rows_with_nan <- sapply(weights_df$Weights, contains_nan)

# Subset the dataframe to exclude rows with NaN in the specified column
weights_df <- weights_df[!rows_with_nan, ]


# Function to list .csv files within a subdirectory
list_csv_files <- function(subdir) {
  list.files(path = file.path(cat.proc, subdir),
             recursive = TRUE,
             pattern = "\\.csv$",
             full.names = FALSE)
}

# List all subfolders within the target density data folder, excluding the target folder itself
subfolders <- list.dirs(path = cat.proc, recursive = FALSE, full.names = FALSE)[-1]

all_geoids <- lapply(subfolders, list_csv_files) 
all_geoids <- unlist(all_geoids)
all_geoids <- as.numeric(gsub("[^0-9]", "", all_geoids))
length(all_geoids)

#-------------------------------------------------------------------------------
# 
# Abb <- c('AL')
# # Use lapply to apply the function to each element of Abb and store the result
# fnms_exist <- lapply(Abb, list_csv_files)
# fnms_exist <- unlist(fnms_exist)
# #length(fnms_exist)
# 
# geoidnms <- unlist(weights_df[1,2]$GEOID_TRACT_10)
# 
# geoidnms = str_pad(geoidnms, 11, pad = "0") 
# fnms = paste0('tract=',geoidnms,'.csv')
# fnms = intersect(fnms, fnms_exist)
# fnms = file.path(cat.proc, Abb, fnms)
# dat1 = readr::read_csv(fnms)
# 
# # average diurnal pattern of the density
# pdat = dat1 %>%
#   mutate(hh = hour(datetime))%>%
#   group_by(hh) %>%
#   summarise(meandensity = mean(agg_density_plpm, na.rm = T))
# 
# maxHH = pdat$hh[pdat$meandensity == max(pdat$meandensity)][1]
# #maxHH # max density hour of the day
# 
# dat = dat1%>%
#   mutate(hh = hour(datetime))%>%
#   group_by(tract_geoid, hh, cnt_tmcs, state)%>%
#   summarise(flow_veh_pspl =  median(agg_flow_plph,na.rm = T)/3600,
#             density_veh_pmpl = median(agg_density_plpm,na.rm = T)/1609.34)%>%
#   ungroup() %>%
#   filter(hh %in% c(7:19)) # this will get us more tracts to work with due to less missing 
# dat = dat %>%
#   dplyr::mutate(GEOID = str_pad(tract_geoid, 11, pad = "0")) %>%
#   dplyr::select(GEOID, hh, density_veh_pmpl)%>%
#   spread(hh, density_veh_pmpl) 
# 
# names(dat)<- c('GEOID',paste0('hh',7:19)) 
# 
# dat$Nna = rowSums(is.na(dat[,-1]))
# 
# # remove the traxts with missing values
# dat = dat %>% dplyr::filter(Nna == 0)
# 

#-------------------------------------------------------------------------------
# Get combined raw density dataframe
load_the_exact_df <- function(file_paths) {
  valid_dfs <- list()  # Create an empty list to store data frames
  
  for (file_path in file_paths) {
    # Read the CSV file
    df <- read.csv(file_path, stringsAsFactors = FALSE)
    
    # Check if the dataframe has exactly 8 columns
    if (ncol(df) == 8) {
      valid_dfs[[length(valid_dfs) + 1]] <- df
    }
  }
  
  # Combine all data frames in the list into a single data frame using dplyr's bind_rows
  combined_df <- bind_rows(valid_dfs)
  
  return(combined_df)
}

# Create a final dataframe to store processed densities
final_df <- data.frame()  # Initialize empty dataframe for final results


# Iteratively process densities for each cbsa
# Assuming 'weights_df' and 'df3' are your dataframes
for (i in 77582:nrow(weights_df)){ 
  #i = 20
  geoidnms <- unlist(weights_df$GEOID_TRACT_10[i])
  Abb <- c()
  for (k in 1:length(geoidnms)){
    Abb <- c(Abb, linker[linker$GEOID == as.numeric(geoidnms[k]), c('st_code')]$st_code)
  }
  
  print(i)
  
  if (length(Abb) == 0) next
  
  if (length(geoidnms) > length(Abb)) {
    # Calculate the difference in length
    length_diff <- length(geoidnms) - length(Abb)
    
    # Repeat the last element of B, length_diff times
    Abb <- c(Abb, rep(Abb[length(Abb)], length_diff))
  }
  
  
  # Abb <- linker %>%
  #   filter(GEOID %in% geoidnms) %>%
  #   pull(st_code) %>%
  #   unique()
  
  fnms_exist <- unlist(lapply(unique(Abb), list_csv_files))
  geoidnms = str_pad(geoidnms, 11, pad = "0") 
  fnms <- paste0('tract=', geoidnms, '.csv')
  tdf <- data.frame(Column1 = fnms, Column2 = Abb)
  fnms <- intersect(fnms, fnms_exist)
  Abb <- tdf[tdf$Column1 %in% fnms, c('Column2')]
  fnms <- file.path(cat.proc, Abb, fnms)
  
  #dat1 = read_csv(fnms)
  dat1 = load_the_exact_df(fnms)
  if (nrow(dat1) == 0) next
  
  # Assuming there are multiple CSVs to read and aggregate
  dat <- lapply(fnms, function(fn) {
    #dat1 <- read_csv(fn)
    dat1 <- dat1 %>%
      mutate(hh = hour(datetime)) %>%
      group_by(tract_geoid, hh, cnt_tmcs, state) %>%
      summarise(flow_veh_pspl = median(agg_flow_plph, na.rm = TRUE) / 3600,
                density_veh_pmpl = median(agg_density_plpm, na.rm = TRUE) / 1609.34,
                .groups = 'drop') %>%
      ungroup() %>%
      filter(hh %in% 7:19) %>%
      mutate(GEOID = str_pad(tract_geoid, 11, pad = "0")) %>%
      select(GEOID, hh, density_veh_pmpl) %>%
      spread(key = hh, value = density_veh_pmpl)
    
    dat1$Nna <- rowSums(is.na(dat1[,-1]))
    dat1 <- dat1 %>% filter(Nna == 0)
    return(dat1)
  })
  
  #print(i)
  
  dat <- as.data.frame(dat[1])
  # if (nrow(dat) == 0) next
 
  dat$Nna <- NULL
  names(dat)<- c('GEOID',paste0('hh',7:19)) 
  
  # For each GEOID, find matching AREALAND_TRACT_10 to calculate weights
  all_matching_arealand <- c()
  for (j in 1:nrow(dat)) {
    #j = 1
    # Extract GEOID for the current row, handling both numeric and character forms
    current_geoid <- as.numeric(dat$GEOID[j])
    
    # Match in 'df3' and compute weights
    matching_arealand <- df3 %>% 
      filter(GEOID_TRACT_10 == current_geoid) %>% 
      pull(AREALAND_TRACT_10)
    
    #matching_arealand <- unique(matching_arealand)
    
    # Skip the iteration if there's no matching row
    if (length(matching_arealand) == 0){
      matching_arealand <- 0
    }
    
    all_matching_arealand <- c(all_matching_arealand, matching_arealand)
  }
  
  weights <- all_matching_arealand / sum(all_matching_arealand)
    
  # Compute weighted sums for the hourly columns
  hourly_columns <- dat[, 2:14]  # Assuming columns 2 to 15 correspond to 7hh to 19hh
  weighted_sum <- sapply(hourly_columns, function(column) sum(column * weights, na.rm = TRUE))
    
  # Create a one-row dataframe for the current row's results
  current_row_df <- data.frame(GEOID = weights_df$GEOID_TRACT_20[i], t(weighted_sum))
  names(current_row_df)[-1] <- paste0((7:19), "hh")
    
  # Append to the final dataframe
  final_df <- rbind(final_df, current_row_df)
  
  
  if (i %% 1000 == 0){
    write.csv(final_df, file = file.path(rdatdir,'new_density2020.csv'))
  }
  
}
  

