# Rscripts to test the spatial coverage by microtypes in 7 cities.

# set working directory to GitHub repo
mywd <- "C:/Users/xiaodanxu/Documents/GitHub/GEMS-data/mfd/partitioning" 
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
auxdatdir <- 'C:/data/CATTLab_delivery/AuxiliaryData' # the blob-microtype ids
#density20 <- 'G:/My Drive/MFD/Tract_level_clustering/National/proc_data/density'
rdatdir <-'C:/data/CATTLab_delivery/National/proc_data/partition_input' # where you put your final results

# if on Ling's Mac: this is where to find the flow density data
# otherwise, set it to your directory where you downloaded the data
cat.proc <- 'C:/data/CATTLab_delivery/National/proc_data/tract_mr16'

# df <- read_csv("D:/Projects/Berkeley/GEMS/typology/us_xwalk_tract_2017_withID_cbsa_only.csv")
# df2 <- read_csv("D:/Projects/Berkeley/GEMS/typology/census_tract_crosswalk_2010_2020.csv")
# df3 <- df2[ ,c('GEOID_TRACT_10', 'AREALAND_TRACT_10')]
# df3 <- unique(df3)


# extracted_geoids <- df[df$cbsa != 99999, 'GEOID']
#subset_df2 <- df2[df2$GEOID_TRACT_20 %in% extracted_geoids$GEOID, c('GEOID_TRACT_20', 'GEOID_TRACT_10', 'AREALAND_TRACT_10')]

# load microtype definition by tract
linker = fread(file.path(auxdatdir,'us_xwalk_tract_2017_withID_cbsa_only.csv')) %>%
  dplyr::select(GEOID, st_code)%>%
  distinct()

# # get weights for matching GEOID2010
# weights_df <- subset_df2 %>%
#   group_by(GEOID_TRACT_20) %>%
#   mutate(Total_AREALAND = sum(AREALAND_TRACT_10)) %>%
#   ungroup() %>%
#   group_by(GEOID_TRACT_20, GEOID_TRACT_10) %>%
#   summarize(Weight = first(AREALAND_TRACT_10) / first(Total_AREALAND), .groups = 'drop') %>%
#   group_by(GEOID_TRACT_20) %>%
#   summarise(GEOID_TRACT_10 = list(GEOID_TRACT_10),
#             Weights = list(Weight))


# contains_nan <- function(x) {
#   if (is.numeric(x)) {
#     return(any(is.nan(x)))  # Returns TRUE if any value is NaN
#   } else {
#     return(FALSE)  # Returns FALSE for non-numeric data (adjust based on your data characteristics)
#   }
# }
# 
# # Apply the function to each element of the list column to identify rows with NaN
# rows_with_nan <- sapply(weights_df$Weights, contains_nan)
# 
# # Subset the dataframe to exclude rows with NaN in the specified column
# weights_df <- weights_df[!rows_with_nan, ]


# Function to list .csv files within a subdirectory
list_csv_files <- function(subdir) {
  list.files(path = file.path(cat.proc, subdir),
             recursive = TRUE,
             pattern = "\\.csv$",
             full.names = FALSE)
}

# List all subfolders within the target density data folder, excluding the target folder itself
subfolders <- list.dirs(path = cat.proc, recursive = FALSE, full.names = FALSE)

all_geoids <- lapply(subfolders, list_csv_files) 
all_geoids <- unlist(all_geoids)
all_geoids <- as.numeric(gsub("[^0-9]", "", all_geoids))
length(all_geoids)

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
for (i in 1:nrow(linker)){ 
  #i = 20
  geoidnms <- unlist(linker$GEOID[i])
  Abb <- c()
  for (k in 1:length(geoidnms)){
    Abb <- c(Abb, linker[linker$GEOID == as.numeric(geoidnms[k]), c('st_code')]$st_code)
  }
  
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
  
  print(i)
  
  dat <- as.data.frame(dat[1])
  # if (nrow(dat) == 0) next
 
  dat$Nna <- NULL
  names(dat)<- c('GEOID',paste0('hh',7:19)) 
  
  # Append to the final dataframe
  final_df <- rbind(final_df, dat)
  
  if (i %% 1000 == 0){
    write.csv(final_df, file = file.path(rdatdir,'density_estimation_2010.csv'),row.names = F)
  }
  
}
  

