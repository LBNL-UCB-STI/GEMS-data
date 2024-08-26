## This code has two sections: Section 1 should be completed before starting the section 2

## Set working directory
mywd <- "D:/Projects/Berkeley/GEMS/" 
setwd(mywd)

#-------------------------------------------------------------------------------
## Install packages
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
p_load(terrainr)
p_load(gridExtra) # plot multiple plots
#p_load(spdep)
remotes::install_version("spdep", version="1.3-4") #SKATER algorithm is prone to produce an error with latest versions.
p_load(clusterSim) 


source('D:/Projects/Berkeley/GEMS/countrywide_partitioning/functions.R')
source('D:/Projects/Berkeley/GEMS/countrywide_partitioning/spatial.partition_functions.R')

#-------------------------------------------------------------------------------
## Section 1: Regionalization and Partitioning
#-------------------------------------------------------------------------------
## Set directories
auxdatdir <- "D:/Projects/Berkeley/GEMS/typology"
rdatdir <-"D:/Projects/Berkeley/GEMS/partition_results_nationwide_2010" ## Where you put your partition results

## User defined parameters: thresholds to stop further partitioning
sub.Ntracts.thr = 2 # if no more than this number of tracts
sub.aland_km2.thr = 4 # or if no more than this number of km2

## Load densities
density_est <- read.csv(file.path(auxdatdir, "density_estimation_2010.csv"))
density_est <- na.omit(density_est) 

## Optional col names
# density_est <- density_est %>%
#   rename(
#     hh7 = X7hh,
#     hh8 = X8hh,
#     hh9 = X9hh,
#     hh10 = X10hh,
#     hh11 = X11hh,
#     hh12 = X12hh,
#     hh13 = X13hh,
#     hh14 = X14hh,
#     hh15 = X15hh,
#     hh16 = X16hh,
#     hh17 = X17hh,
#     hh18 = X18hh,
#     hh19 = X19hh
#   )

## Load cbsa data. Use the correct CBSA tract information: 2010 or 2020
df <- read_csv(file.path(auxdatdir, 'us_xwalk_tract_2017_withID_cbsa_only.csv'))
cbsa_df <- unique(df[df$cbsa != 99999, 'cbsaname'])

## Iterate over CBSAs
for(y in 1:length(cbsa_df$cbsaname)){
  if (y %in% c(480, 713)) {next}
  
  print(y)
  
  city.sel0 = cbsa_df$cbsaname[y]
  city.sel = gsub("/", "-",   city.sel0)
  
  words <- unlist(strsplit(city.sel, ", "))
  city = words[1]
  state = words[2]
  
  city_tracts = df %>%
    #left_join(citymap) %>%
    drop_na() %>%
    filter(cbsaname != "") %>%
    filter(cbsaname %in% city.sel0)
  
  geoidnms = city_tracts$GEOID
  geoidnms = str_pad(geoidnms, 11, pad = "0") # the filename has 11 digits for tractid
  state.sel = substr(geoidnms[1],1,2)
  tracts = tigris::tracts(state = state.sel, year = 2010) %>% 
    dplyr::filter(as.integer64(GEOID10) %in% as.integer64(city_tracts$GEOID[city_tracts$cbsaname ==  city.sel0])) %>%
    mutate(aland_km2 = ALAND10/1000000) %>%
    filter(aland_km2<10) %>%
    mutate(tract_geoid = as.integer64(GEOID10))# only keep the tracts less than 10km2
  
  tracts$GEOID10 <- as.double(tracts$GEOID10)
  names(tracts)[names(tracts) == "GEOID10"] <- "GEOID"
  
  tracts1 <- left_join(tracts, density_est, by = "GEOID")
  tracts1 <- tracts1 %>% drop_na()
  
  ## Regions with less than or equal to 2 tracts are not clustered..
  if (nrow(tracts1) <= 2){
    data_path = file.path(rdatdir, paste0(state))
    ## Check if the folder exists
    if (!dir.exists(data_path)) {
      ## Create the folder
      dir.create(data_path, recursive = TRUE)  ## recursive=TRUE allows creating nested directories
    } 
    save(tracts1, file = file.path(paste0(data_path,'/', city,'.partition.results.RData')))
    next
  }
  
  ## Max density hour of the day
  columns_of_interest <- paste0("hh", 7:19)
  column_sums <- colSums(st_drop_geometry(tracts1[columns_of_interest]), na.rm = TRUE)
  
  ## Identify the column name with the maximum sum
  maxcol <- names(column_sums)[which.max(column_sums)]
  maxHH <- as.numeric(gsub("\\D", "", maxcol))
  
  ## Now partition the connected region into homogeneous zones
  # Need to find connected tracts
  poly1 <- sf::st_union(tracts1) %>% # dissolve boundaries between polygons
    sf::st_cast(to = 'POLYGON') %>% # recast multi-polygon sf object into sf object of multiple single polygons
    sf::st_sf() %>% # turn into sf object for easier manipulation
    mutate(regID = row.names(.)) # create unique ID for each merged blob
  unique(poly1$regID)
  
  kk = st_intersection(tracts1, poly1) %>%
    dplyr::select(GEOID, regID) %>%
    st_drop_geometry() %>%
    distinct()
  
  ## Add the region definition to the original data
  tracts1 = tracts1 %>%
    left_join(kk) %>%
    group_by(regID) %>% ## Compute number of tracts, and land area in each region
    mutate(reg_aland_km2 = sum(aland_km2,na.rm = T),
           reg_Ntracts = n())%>%
    ungroup() %>%
    mutate(cluster = '00') ## Initialize the cluster labeling
  
  ## Define thresholds of number of tracts and land area to warrant initial partition
  Ntracts.thr = 4
  aland_km2.thr = 1 
  ## if number of tracts <= Ntract.thr or land area <= aland_km2.thr,  satisfied, no partition for the region
  regs.sel = unique(tracts1$regID[tracts1$reg_aland_km2 > aland_km2.thr & tracts1$reg_Ntracts > Ntracts.thr])
  
  ## First loop over the eligible regions to conduct first round of partition based on DBI-determined k and skater
  for(regid in regs.sel){
    #regid = regs.sel[7]
    init_partition(regid) ## Perform initial partitioning for each continous region
  }
  
  ## Obtain simple features
  tracts1 <- st_as_sf(tracts1)
  
  ## Within each cluster by a regID, initialize a subcluster id
  tracts1$sub.cluster = tracts1$cluster
  
  ## Now compute the land area and number of tracts of each partition
  tracts1 <- tracts1 %>%
    group_by(regID, sub.cluster)%>%
    mutate(cl.aland_km2 = sum(aland_km2,na.rm = T),
           cl.Ntracts = n())%>%
    ungroup()
  
  tracts1 = unite(tracts1, col = 'unitedID', c('regID','sub.cluster'), sep='-' , remove = F) 
  
  more.part = tracts1%>%
    filter(cl.aland_km2 > sub.aland_km2.thr & cl.Ntracts > sub.Ntracts.thr)%>%
    dplyr::select(unitedID)%>%
    st_drop_geometry()%>%
    distinct()
  N.more = dim(more.part)[1]
  
  n = 0
  while(N.more > 0 & n < 10000){
    
    for(idx in more.part$unitedID){
      if ((y == 68 && idx %in% c("19-01", "bad clusters")) || 
          (y == 339 && idx %in% c("5-02", "bad clusters")) ||
          (y == 375 && idx %in% c("5-02", "4-02", "12-02", "bad clusters")) ||
          (y == 618 && idx %in% c("15-02", "15-04", "bad clusters"))) {
        next # Skip the bad clusters
      }
      print(idx)
      split2_partition(idx) ## Split the partiion into 2 using skater
    }
    
    ## Recompute the area and ntracts of each partiion
    tracts1 <- tracts1 %>%
      group_by(regID, sub.cluster)%>%
      mutate(cl.aland_km2 = sum(aland_km2,na.rm = T),
             cl.Ntracts = n())%>%
      ungroup()
    
    tracts1 = unite(tracts1, col = 'unitedID', c('regID','sub.cluster'), sep='-' , remove = F) 
    
    more.part = tracts1%>%
      filter(cl.aland_km2 > sub.aland_km2.thr & cl.Ntracts > sub.Ntracts.thr)%>%
      dplyr::select(unitedID)%>%
      st_drop_geometry()%>%
      distinct()
    N.more = dim(more.part)[1]
    n = n + 1
  }
  
  data_path = file.path(rdatdir, paste0(state))
  # Check if the folder exists
  if (!dir.exists(data_path)) {
    # Create the folder
    dir.create(data_path, recursive = TRUE)  ## recursive=TRUE allows creating nested directories
  }  
  
  save(tracts1, file = file.path(paste0(data_path,'/', city,'.partition.results.RData')))
  
}



#-------------------------------------------------------------------------------
## Section 2: Merging orphan tracts 
#-------------------------------------------------------------------------------

# List all .RData files recursively
rdata_files <- list.files(path = rdatdir, pattern = "\\.RData$", full.names = TRUE, recursive = TRUE)

# Create dataframe and matrix
df_results <- as.data.frame(matrix(NA, nrow = 930, ncol = 3))

column_names <- c("cbsa", "n_orphans", "dbi_increase")
colnames(df_results) <- column_names

final_results <- as.data.frame(matrix(NA, nrow = 1, ncol = 3))
colnames(final_results) <- c('GEOID', 'unitedID', 'cbsa')

k = 0
## Iterate over CBSAs
for(file in rdata_files){
 
  #file <- rdata_files[309]
  load(file)
  
  ## Extract the cbsa
  if (k == 310){
    cbsa_id = 31140 #"Louisville-Jefferson County, KY-IN"
  }

  if (k != 310){
    city <- sub(".*/([^/]+)/([^/]+)\\.partition\\.results\\.RData", "\\2, \\1", file)
    cbsa_id <- unique(subset(df, cbsaname == city)$cbsa)
    df_results[k,1] <- city
  }

  if(all(c("cl.aland_km2", "cl.Ntracts") %in% names(tracts1)) == TRUE & nrow(tracts1) >= 3){
    
    tracts1 <- tracts1[!duplicated(tracts1$GEOID), ]
    ## Find all the neighbors
    #neighbors <- st_touches(tracts1, tracts1, sparse = T)
    ## Get the whole subset of orphan tracts
    s.partitions.pool <- subset(tracts1, cl.Ntracts %in% 1)
    
    ## Find the distribution of the land area of the orphan tracts
    #s.area <- s.partitions.pool$cl.aland_km2
    partition.sizes <- c(3.0)
    #partition.sizes <- seq(0.5, 5.0, by = 0.5)
    #partition.sizes <- partition.sizes[partition.sizes <= max(s.area)]
    
    ## Create an empty matrix to store result
    output <- matrix(ncol=3, nrow=1)
    
    ## Hyper-parameter tuning using all percentiles of the area of orphans
    for (iter in 1:length(partition.sizes)){
      ## Assign initial L2 distance between an orphan and a candidate
      #v <- 0.00001
      ## Find the land area threshold for this iteration
      area.thd <- partition.sizes[iter]
      tracts2 <- tracts1
      
      ## store the updated unitedID
      tracts2$unitedID.new <- tracts1$unitedID
      ## Get a subset of orphan tracts satisfying the above land area threshold
      s.partitions <- s.partitions.pool[s.partitions.pool$cl.aland_km2 < area.thd,] 
      
      print(nrow(s.partitions))
      if (nrow(s.partitions) == 0){
        next
      }
      
      # Merge orphans for a given area percentile 
      results <- merge_orphans(tracts2, area.thd)
      # n <- as.numeric(results$n_orphans)
      # tracts2 <- data.frame(results$updated_tracts)
      
      n <- as.numeric(results[1])
      # the resulting tracts2
      tracts2 <- data.frame(results[2])
      dft <- data.frame(results[3])
      #print(n) 
      
      output[iter,1] <- area.thd
      output[iter,2] <- n
      output[iter,3] <- (Davies_Bouldin_2(tracts2) - Davies_Bouldin_1(tracts1))/Davies_Bouldin_1(tracts1)
      colnames(output) <- c('Area', 'N.orphans', 'DBI')
    }
    
    # Insert matrix columns into dataframe
    df_results[k, 2] <-  output[1, 2]
    df_results[k, 3] <-  output[1, 3]
    # df_results[k, 2:11] <-  output[, 2]
    # df_results[k, 12:21] <-  output[, 3]
    
    # Update `unitedID.new` values to NA where partition did not start
    tracts2$unitedID.new[grepl("^\\d{1,2}-\\d{2}$", tracts2$unitedID.new)] <- NA
    
    # Update `unitedID.new` values to NA where cl.aland_km2 >= 4 and cl.Ntracts >= 3
    tracts2$unitedID.new[tracts2$cl.aland_km2 >= 4 & tracts2$cl.Ntracts > 2] <- NA
    
    subdf <- tracts2[c('GEOID', 'unitedID.new')]
    subdf <- st_drop_geometry(subdf)
    cbsa_vec <- rep(cbsa_id, nrow(tracts2)) 
    
  }
  
  if (all(c("cl.aland_km2", "cl.Ntracts") %in% names(tracts1)) == TRUE & nrow(tracts1) < 3){
    
    tracts1$unitedID[grepl("^\\d{1,2}-\\d{2}$", tracts1$unitedID)] <- NA
    tracts1$unitedID[tracts1$cl.aland_km2 >= 4 & tracts1$cl.Ntracts > 2] <- NA
    
    subdf <- tracts1[c('GEOID', 'unitedID')]
    subdf <- st_drop_geometry(subdf)
    cbsa_vec <- rep(cbsa_id, nrow(tracts1)) 
    
  }
  
  if (all(c("cl.aland_km2", "cl.Ntracts") %in% names(tracts1)) == FALSE){
    
    subdf <- tracts1[c('GEOID')]
    subdf <- st_drop_geometry(subdf)
    uid <- rep('NA', nrow(tracts1))
    subdf$unitedID <- uid
    cbsa_vec <- rep(cbsa_id, nrow(tracts1)) 
    
  }
  
  subdf$cbsa <- cbsa_vec
  colnames(subdf) <- c('GEOID', 'unitedID', 'cbsa')
  
  if (nrow(subdf) > 0){
    final_results <- rbind(final_results, subdf)
  }
  
  print(paste("cbsa:", k))
  k <- k + 1
  
}

## Save final results
write.csv(df_results, "D:/Projects/Berkeley/GEMS/final_merginging_orphans_results.csv", row.names = FALSE)

final_results <- final_results[-1,]
write.csv(final_results, "D:/Projects/Berkeley/GEMS/final_partition_results.csv", row.names = FALSE)


#-------------------------------------------------------------------------------

