# LJ 8/10/2024
# Rscript to process the partition results from Duleep, 
# assigning unique partition ids and generate cross walk, 
# only keep the partitions containing >=2 tracts, which will need re-aggregation of flow and density.
# link to repo
mywd = "C:/Users/xiaodanxu/Documents/GitHub/GEMS-data/mfd"
#mywd <- "/Users/lingjin/Dropbox/Research/FHWA_GeoType/Work/Rscripts"
setwd(mywd)
#source('initialization.R')
#source('functions.R')
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



#source('initialization.R')
source('functions.R')

# change to data dir
#mywd = "C:/data"
#mywd <- "/Users/lingjin/Dropbox/Research/FHWA_GeoType/Work/Rscripts"
#setwd(mywd)
auxdatdir <- 'C:/data/CATTLab_delivery/AuxiliaryData' # the old microtypes
rdatadir <-  "C:/data/RData/MFD_v2"
figuredir <- "C:/data/Figures/MFD_v2"
#outdir  <-  "C:/data/"

cat.downloads <- 'C:/data/CATTLab_delivery/National/'


part.dir <- 'C:/data/CATTLab_delivery/National/proc_data/partition_input'
# part.dir <- 'C:/RData/MFD_v2/partition_input'

linker = fread(file.path(auxdatdir,'ccst_geoid_key_tranps_geo_with_imputation.csv'))%>% # old mcirotype xwalk
  select(GEOID, microtype, st_code, cbsa, cbsaname,geotype)

##################

pt <- fread(file.path(part.dir,'final_partition_results_2010.csv')) %>%
  filter(!is.na(unitedID)) %>%
  unite('partID',unitedID:cbsa,sep = "_", remove = F)

# get some summaries:

sum.pt <- pt %>%
  group_by(partID)%>%
  summarise(ntracts = n())

table(sum.pt$ntracts)
round(prop.table(table(sum.pt$ntracts)),2)

# remove single tract partitions, b/c they do not need further reaggregation
pt <- pt %>%
  group_by(partID) %>%
  mutate(ntracts = n()) %>%
  filter(ntracts >1) %>%
  select(-ntracts)

pt <- pt %>%
  left_join(linker)

# create readable partition id
pid <- pt %>%
  select(partID) %>%
  distinct()
pid$partid <- 1:(dim(pid)[1])

pt <- pt %>%
  left_join(pid)


# save the tracts and xwalk that need re-aggregation.

write.csv(pt, file = file.path(part.dir,paste0('reagg_partition_xwalk_092424.csv')),row.names = F)
save(pt, file = file.path(part.dir, 'reagg_partition_xwalk_092424.RData'))

pt.agg <- pt %>%
  mutate(tract_geoid = as.numeric(GEOID)) %>%
  select(tract_geoid, cbsa, partid, cbsaname, st_code)

# save all the partition xwalk including the ones with single tract
pt1 <- fread(file.path(part.dir,'final_partition_results_2010.csv')) %>%
  filter(!is.na(unitedID)) %>%
  unite('partID',unitedID:cbsa,sep = "_", remove = F) %>%
  group_by(partID) %>%
  mutate(ntracts = n()) %>%
  ungroup() %>%
  left_join(linker) %>%
  filter(ntracts == 1) %>%
  select(-ntracts)%>%
  mutate(tract_geoid = as.numeric(GEOID)) %>%
  select(tract_geoid, cbsa, cbsaname, st_code)



pt1$partid = (1:(dim(pt1)[1]))+max(unique(pt.agg$partid)) # assign these single tract partitions to unique partition ids. 

save(pt1, file = file.path(part.dir, 'single_tract_partition_xwalk_092424.RData'))

pt = bind_rows(pt1, pt.agg)

save(pt, file = file.path(part.dir, 'all_partition_xwalk_092424.RData'))
