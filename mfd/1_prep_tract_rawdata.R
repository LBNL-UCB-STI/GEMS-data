# LJ 8/1/2024 Rscripts for MFD on new partitions
# prepare tract level raw here data by replacing the microtype 3 data aggregating 1-5frc with Kaveh's new data aggregating 3-5frc
# redo aggregation to the partition level
# extract the upper bound for next MFD parameter esitmation

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

auxdatdir <- 'C:/data/CATTLab_delivery/AuxiliaryData' # the blob-microtype ids
#rdatadir <-  "../RData/CATTLab/National_Batch/ApproachA"
#tabdir <- '../Tables/CATTLab/National_Batch/ApproachA'
#figuredir <- "../Figures/CATTLab/National_Batch/ApproachA"
#outdir  <-  "../Data/"
part.dir <- 'C:/data/CATTLab_delivery/National/proc_data/partition_input' # partition xwalk


cat.downloads <- 'C:/data/CATTLab_delivery/National/'

cat.proc = 'C:/data/CATTLab_delivery/National/proc_data/all_tracts_frc3_5'  # tract level data aggregating frc3-5 consistently
cat.proc2 <- 'C:/data/CATTLab_delivery/National/proc_data/m3_tracts_frc1_5' # only tract level data for microtype 3 aggregating frc1-5

batch.dirs <- c('OneDrive_2_2-1-2023','OneDrive_3_2-1-2023','OneDrive_4_2-1-2023','OneDrive_5_2-1-2023','OneDrive_6_2-1-2023') # tracts with microtype 3 aggregating frc1-5

batch.dirs2 <- file.path(cat.downloads,'Microtype3Tracts_FRC345_8-31-2023/20230828/agg_microtype3_frc3plus') # microtype 3 aggregating frc3-5

approach = 'ApproachA'

# get the lanemiles ready for each tract
lms = fread(file = file.path(auxdatdir,'combined_HERE_lanemiles.csv')) %>%
  mutate(frc345_lanemiles = frc3_lanemiles + frc4_lanemiles + frc5_lanemiles) %>%
  select(tract_geoid,frc345_lanemiles)%>%
  mutate(tract_geoid = as.numeric(tract_geoid))

#1. get tract-microtype xwalk

linker = fread(file.path(auxdatdir,'ccst_geoid_key_tranps_geo_with_imputation.csv'))
linker.m3 = as.numeric(unique(linker$GEOID[linker$microtype == 3]))

landarea = fread(file.path(auxdatdir,'row.csv')) %>%
  mutate(aland_km2 = aland/10^6) %>%
  rename(GEOID = tract)%>%
  select(GEOID, aland_km2)
# xwalk
# xwalk file, and pad 0 in GEOID
xwalk = fread(file.path(auxdatdir, 'us_xwalk_tract_2017.csv')) %>%
  rename(GEOID = tract) %>% # %>% # an integer64 class
  #  mutate(GEOID = as.character(str_pad(GEOID, width=11, side="left", pad="0")))
  #  filter(state == 'California') %>% # get california only
  left_join(landarea) %>%
  left_join(linker) %>% # FID is the blob ID
  select(FID, st_code, cty, ctyname, GEOID, cbsa, cbsaname, aland_km2, microtype, geotype, spatial_id)%>%
  rename(blob_id = FID) 

xwalk.tract = xwalk %>%
  select(blob_id,GEOID, st_code,cbsa,microtype, cbsaname, geotype, aland_km2)%>%
  distinct()%>%
  rename(tract_geoid = GEOID)

# 2. extract by tract the new microtype 3 data that aggregated frc3-5

# function to extract data
tract_extractm3 <- function(fnm, st){
  library(pacman)
  p_load(R.utils) # for reading gz files
  p_load(readr) # for read files from a zip foler
  
  #get the tract id
  tid = substr(fnm, nchar("tract=")+1, nchar(fnms)-nchar(".csv"))
  dat <- read_csv(file.path(m3dir,fnm)) %>% 
    select(datetime, cnt_tmcs,contains("V1") )%>%
    rename(agg_flow_plph = V1_aggflow_plph,
           agg_density_plpm = V1_aggdensity_plpm,
           aggspeed_mph = V1_aggspeed_mph )%>%
    mutate(tract_geoid = as.numeric(tid),
           state = st)
  
  if(length(dat)>0){
    return(dat)
  }else{return(NULL)}
  
}

m3dirs = list.dirs(batch.dirs2, recursive = F)

# loop over state and extract tract data

for(m3dir in m3dirs){
  state <- substr(m3dir, nchar(m3dir)-1, nchar(m3dir))
  fnms = list.files(m3dir,recursive = F)
  
  registerDoParallel(cores = 15)
  tic()
  res <- foreach(k=1:length(fnms),
                     .combine=rbind,
                     .packages = c('dplyr','tidyr','vroom','data.table','pacman','R.utils','readr')
  )  %dopar% {
    tract_extractm3(fnms[k], st=state)
  }
  toc()
  
  
  stopImplicitCluster()
  
  write.csv(res, file = file.path(cat.proc,paste0('m3_tract_frc345_state_',state,'.csv')),row.names = F)
  
}


# 2.2. extract 20 percentile upper bounds from the updated microtype 3 tracts

m3resdirs = list.files(cat.proc, recursive = F)
upperbound_pct = 0.2  # the top 30 percentile as upper bound


for(m3res in m3resdirs){
  master.dat = fread(file = file.path(cat.proc, m3res))
  master.dat$microtype = 3
  
  cdat = master.dat %>%
    group_by(tract_geoid, state, microtype) %>% # error handling and filter out tracts with bad or limited data
    mutate(nsample = n(),
           mdensity = median(agg_density_plpm , na.rm = T),
           mflow = median(agg_flow_plph, na.rm =T)) %>%
    filter(nsample >=1000,
           mdensity >0.001,
           mflow > 0.001) %>%
    select(-nsample, -mdensity, -mflow)%>%
    mutate(maxdens = max(agg_density_plpm,na.rm = T),
           mindens = min(agg_density_plpm, na.rm = T))%>%
    mutate( dens.bins = cut(agg_density_plpm, 
                            breaks = seq(unique(mindens),unique(maxdens),
                                         by = (unique(maxdens)-unique(mindens))/60))) %>%
    ungroup()%>%
    group_by(tract_geoid,microtype, state, dens.bins) %>%
    mutate(c.flow.l = quantile(agg_flow_plph, (1-upperbound_pct)),
           c.flow.h = quantile(agg_flow_plph,0.99),
           nsample = n()) %>%
    ungroup() %>%
    filter((agg_flow_plph >= c.flow.l & agg_flow_plph <= c.flow.h ) | nsample <=3 ) %>%
    dplyr::select(-c.flow.l, -c.flow.h, -nsample, -maxdens, -mindens, -dens.bins)
  
write.csv(cdat, file = file.path(cat.proc,paste0(substr(m3res, 1, str_length(m3res)-4), 'upperbounds.csv')),row.names = F)
}

# 3. extract by tract the other microtype data that aggregated frc3-5, replacing microtype 3 tracts with the updated data from step 2 above
#     also extract upper 20 percentile of each tract
# save the data for each state for later MFD estimates


for(j in 1:length(batch.dirs)){ # loop over the zipped folders that contains multiple states
  
  
  batch.dir = file.path(cat.downloads, batch.dirs[j])
  catdatdirs = list.dirs(batch.dir, recursive = F)
  states = substr(catdatdirs,nchar(catdatdirs)[1]-1, nchar(catdatdirs)[1])
  
  
  for(i in 1:length(states)){ # loop over states
    #  state = 'CA' # state
    spunit = 'tract' # spatial unit
    approach = 'ApproachA' # approach A or C
    state = states[i]
    
    print(state)
    
    # load the corrected microtype 3 upper bound data for this state
    
    m3dat = fread (file = file.path(cat.proc,paste0('m3_tract_frc345_state_',state,'upperbounds.csv')))
    
    # dirname = file.path(datdir, state, paste0(state,"_",spunit,"_aggfd_data"),approach) # flow density data
    # 
    #   fnms = file.path(dirname,as.vector(list.files(dirname)))
    
    testdir = file.path(catdatdirs[i],paste(states[i],spunit,'aggfd_data.zip',sep='_'))
    fnms = unzip(zipfile = testdir, list = TRUE)
    fnms = fnms %>% filter(grepl(approach,Name), Length>0)
    
    #  df1 <- read_csv(unzip(testdir,fnms$Name[1]))
    
    
    registerDoParallel(cores = 15)
    tic()
    res <- foreach(k=1:length(fnms$Name),
                       .combine=rbind,
                       .packages = c('dplyr','tidyr','vroom','data.table','pacman','R.utils','readr')
    )  %dopar% {
      print(k)
      # only extract microtype 12
      # extract only if tract is not microtype 3
     mfd.proc.extractv1_exclude3( spatial_unit = spunit, upperbound_pct = 0.2, zipTF = T, zipfolder=testdir, fnm.inzip=fnms$Name[k])
    } # end dopar
    
    toc()
    
    
    stopImplicitCluster()
    
    if(length(res)>0){
      res$state = state
      res$tract_geoid = as.numeric(res$tract_geoid)
      m3dat$tract_geoid = as.numeric(m3dat$tract_geoid)
      res = bind_rows(res, m3dat)
      write.csv(res, file = file.path(cat.proc, paste0('bytract_upperbounds_state_',state, '.csv')),row.names = F)
    }
        
    
    
  } # end loop state
  
} # end loop zip folders





# 4. extract the tracts that appeared in the partition of 2 or more tracts, these tracts will need to retain all the data so that they can be re-aggregated.

# load the partition data
pt <- fread(file = file.path(part.dir,paste0('reagg_partition_xwalk.csv'))) %>%
  mutate(tract_geoid = as.numeric(GEOID)) %>%
  select(tract_geoid, cbsa, cbsaname, partid, st_code)





# loop over states to accumulate the tracts that goes into partition and get upper bounds

upperbound_pct = 0.2 
res.all = NULL # all data
res.upperbound = NULL # only save the upper bounds

for(j in 1:length(batch.dirs)){ # loop over the zipped folders that contains multiple states
  
  
  batch.dir = file.path(cat.downloads, batch.dirs[j])
  catdatdirs = list.dirs(batch.dir, recursive = F)
  states = substr(catdatdirs,nchar(catdatdirs)[1]-1, nchar(catdatdirs)[1])
  
  
  for(i in 1:length(states)){ # loop over states
    #  state = 'CA' # state
    spunit = 'tract' # spatial unit
    approach = 'ApproachA' # approach A or C
    state = states[i]
    
    print(state)
    
    # load the corrected microtype 3 all data for this state
    
    m3dat = fread (file = file.path(cat.proc,paste0('m3_tract_frc345_state_',state,'.csv')))%>%
      mutate(tract_geoid = as.numeric(tract_geoid))%>%
      filter(tract_geoid %in% pt$tract_geoid)%>%
      mutate(microtype = 3)
    
    # also accumulate the tract data from other microtypes
    testdir = file.path(catdatdirs[i],paste(states[i],spunit,'aggfd_data.zip',sep='_'))
    fnms = unzip(zipfile = testdir, list = TRUE)
    fnms = fnms %>% filter(grepl(approach,Name), Length>0)
    
    #  df1 <- read_csv(unzip(testdir,fnms$Name[1]))
    
    
    registerDoParallel(cores = 15)
    tic()
    res <- foreach(k=1:length(fnms$Name),
                   .combine=rbind,
                   .packages = c('dplyr','tidyr','vroom','data.table','pacman','R.utils','readr')
    )  %dopar% {
      print(k)
         # extract only if tract is not microtype 3
   proc.extractv1_exclude3( selecttracts = pt$tract_geoid[pt$st_code == state], zipTF = T, zipfolder=testdir, fnm.inzip=fnms$Name[k])
    } # end dopar
    
    toc()
    
    
    stopImplicitCluster()
    
    if(length(res)>0){
      res$state = state
      res = bind_rows(res, m3dat)
    
    # re-aggregate according to the partitions
    res <- res %>%
      left_join(pt)%>%
      left_join(lms) # merge in lane-miles
    
    # redo aggregation
    resdat = res %>%
      group_by(partid,cbsa,cbsaname,state,datetime)%>%
      summarise(flow_plph = sum(agg_flow_plph * frc345_lanemiles)/sum(frc345_lanemiles),
                density_plpm = sum(agg_density_plpm * frc345_lanemiles)/sum(frc345_lanemiles) )%>%
      rename(agg_flow_plph = flow_plph,
             agg_density_plpm = density_plpm )%>%
      mutate(aggspeed_mph = agg_flow_plph/agg_density_plpm)
    
    # get the upper bounds
    cdat = resdat %>%
      group_by(partid, cbsa, cbsaname, state) %>%
      mutate(nsample = n()) %>%
      filter(nsample >=1000) %>%
      select(-nsample) # only process if tracts have more than 1000 points
    if(nrow(cdat)>0){
    cdat = cdat %>%
      group_by(partid, cbsa, cbsaname, state)%>%
      mutate(maxdens = max(agg_density_plpm,na.rm = T),
             mindens = min(agg_density_plpm, na.rm = T))%>%
      mutate( dens.bins = cut(agg_density_plpm, 
                              breaks = seq(unique(mindens),unique(maxdens),
                                           by = (unique(maxdens)-unique(mindens))/60))) %>%
      ungroup()%>%
      group_by(partid, cbsa, cbsaname, state, dens.bins ) %>%
      mutate(c.flow.l = quantile(agg_flow_plph, (1-upperbound_pct)),
             c.flow.h = quantile(agg_flow_plph,0.99),
             nsample = n()) %>%
      ungroup() %>%
      filter((agg_flow_plph >= c.flow.l & agg_flow_plph <= c.flow.h ) | nsample <=3 ) %>%
      dplyr::select(-c.flow.l, -c.flow.h, -nsample, -maxdens, -mindens, -dens.bins)
   
     res.upperbound = rbind(res.upperbound, cdat)
    
    } # only process if cdat nonempty
    
    # accumulate the final results
    res.all = rbind(res.all,resdat)
    
    
    
    } # if there is res
    
    
  } # end loop state
  
} # end loop zip folders


 save(res.all, file = file.path(cat.proc, 'reaggreg_partition_alldat.RData'))
 save(res.upperbound, file = file.path(cat.proc, 'reaggreg_partition_upperbound.RData'))
 

 # 5. extract by tract the old microtype 3 that aggregated frc1-5, 
 # 
 # save the m3 data
 
 for(j in 1:length(batch.dirs)){ # loop over the zipped folders that contains multiple states
   
   
   batch.dir = file.path(cat.downloads, batch.dirs[j])
   catdatdirs = list.dirs(batch.dir, recursive = F)
   states = substr(catdatdirs,nchar(catdatdirs)[1]-1, nchar(catdatdirs)[1])
   
   
   for(i in 1:length(states)){ # loop over states
     #  state = 'CA' # state
     spunit = 'tract' # spatial unit
     approach = 'ApproachA' # approach A or C
     state = states[i]
     
     print(state)
     

     # dirname = file.path(datdir, state, paste0(state,"_",spunit,"_aggfd_data"),approach) # flow density data
     # 
     #   fnms = file.path(dirname,as.vector(list.files(dirname)))
     
     testdir = file.path(catdatdirs[i],paste(states[i],spunit,'aggfd_data.zip',sep='_'))
     fnms = unzip(zipfile = testdir, list = TRUE)
     fnms = fnms %>% filter(grepl(approach,Name), Length>0)
     
     #  df1 <- read_csv(unzip(testdir,fnms$Name[1]))
     
     fnms$tract_geoid = as.numeric(substr(fnms$Name,nchar(fnms$Name)-14, nchar(fnms$Name)-4))
     
     fnms = fnms %>%
       filter(tract_geoid %in% linker.m3)  # only select files that are microtype3 tracts
     
     registerDoParallel(cores = 15)
     tic()
     res <- foreach(k=1:length(fnms$Name),
                    .combine=rbind,
                    .packages = c('dplyr','tidyr','vroom','data.table','pacman','R.utils','readr')
     )  %dopar% {
       print(k)
 
       proc.extractv1(zipTF = T, zipfolder=testdir, fnm.inzip=fnms$Name[k])
     } # end dopar
     
     toc()
     
     
     stopImplicitCluster()
     
     if(length(res)>0){
       res$state = state
       res$tract_geoid = as.numeric(res$tract_geoid)
       write.csv(res, file = file.path(cat.proc2, paste0('bytract_m3frc1-5_state_',state, '.csv')),row.names = F)
     }
     
     
     
   } # end loop state
   
 } # end loop zip folders
 


 # 6. extract by tract the upper bound from old microtype 3 that aggregated frc1-5, 
 # 
 # save the m3 data
 
 m3resdirs = list.files(cat.proc2, recursive = F)
 upperbound_pct = 0.2  # the top 20 percentile as upper bound
 
 
 for(m3res in m3resdirs){
   master.dat = fread(file = file.path(cat.proc2, m3res))
   # master.dat$microtype = 3
   
   cdat = master.dat %>%
     group_by(tract_geoid, state, microtype) %>% # error handling and filter out tracts with bad or limited data
     mutate(nsample = n(),
            mdensity = median(agg_density_plpm , na.rm = T),
            mflow = median(agg_flow_plph, na.rm =T)) %>%
     filter(nsample >=1000,
            mdensity >0.001,
            mflow > 0.001) %>%
     select(-nsample, -mdensity, -mflow)%>%
     mutate(maxdens = max(agg_density_plpm,na.rm = T),
            mindens = min(agg_density_plpm, na.rm = T))%>%
     mutate( dens.bins = cut(agg_density_plpm, 
                             breaks = seq(unique(mindens),unique(maxdens),
                                          by = (unique(maxdens)-unique(mindens))/60))) %>%
     ungroup()%>%
     group_by(tract_geoid,microtype, state, dens.bins) %>%
     mutate(c.flow.l = quantile(agg_flow_plph, (1-upperbound_pct)),
            c.flow.h = quantile(agg_flow_plph,0.99),
            nsample = n()) %>%
     ungroup() %>%
     filter((agg_flow_plph >= c.flow.l & agg_flow_plph <= c.flow.h ) | nsample <=3 ) %>%
     dplyr::select(-c.flow.l, -c.flow.h, -nsample, -maxdens, -mindens, -dens.bins)
   
   write.csv(cdat, file = file.path(cat.proc2,paste0(substr(m3res, 1, str_length(m3res)-4), 'upperbounds.csv')),row.names = F)
 }
 