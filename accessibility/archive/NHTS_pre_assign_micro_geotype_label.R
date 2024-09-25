## Rscript to compile travel outcomes from NHTS 2017
# Ling Jin 5/6/2020
# 6/3/2020 adding 2009 data
# merge trip o-d to tract ids and then with micro/geotypes
# 7/8/2020: match to national data
# 3/15/2021: use the transgeo version of clustering results


mywd = "C:/FHWA/For FHWA folks/Mode_choice_input_generation"
setwd(mywd)

 # nhts.dir <- 'C:/NHTS/2017_tripct/'

source('./Code/initialization.R')
source('./Code/functions.R')
# figuredir <- "../Figures/"
# outdir  <-  "../Data/"
# rdatadir <-  "../RData"
datdir = './Data/'

library(tidyverse)
library(reshape2)

# merged transgeo microtypes (remember to merge last b/c there are NA labels)
# transgeof  = 'contiguous_microtype_blobs_key_labeled.csv' ->Xiaodan's notes: cannot find this file, tried the other version instead
transgeof  = 'ccst_geoid_key_transp_geo_with_imputation.csv'
# read in transgeo clustering results
linker = fread(file.path(datdir, transgeof)) %>%
  dplyr::select(GEOID, microtype, geotype) %>%
  mutate(GEOID = as.numeric(GEOID))


### functions #########
match.microtype <- Vectorize(function(geoid){
  tt = as.numeric(linker[linker$GEOID == geoid,'microtype'])
  return(tt)
})

match.geotype <- Vectorize(function(geoid){
  tt = as.numeric(linker[linker$GEOID == geoid,'geotype'])
  return(tt)
})

#match.microtype <- function(geoid){
#  tt = as.numeric(linker[linker$GEOID == geoid,'cluster6'])
#  return(tt)
#}

#match.geotype <- function(geoid){
#  tt = as.numeric(linker[linker$GEOID == geoid,'cluster9'])
#  return(tt)
#}


#### 1. process the trip o-d locations  ######
triploc = fread(file = file.path(datdir,'tripct.csv') ) %>%
  unite('o_geoid',ORIG_ST : ORIG_CT, sep='') %>%
  unite('d_geoid',DEST_ST : DEST_CT, sep = '') %>%
  mutate(o_geoid = as.numeric(o_geoid),
         d_geoid = as.numeric(d_geoid)) 

#### 2. merge to trip data to get trip purposes ######
trips = fread(file = file.path(datdir,'trippub.csv'))
colnames(trips) = tolower(colnames(trips))
names(trips)
dim(trips) # 923572 trips

#rename to CA convention
trips = trips %>%
  rename(sampno = houseid,
         perno = personid)
trips = trips %>% 
  filter(!is.na(sampno) & !is.na(perno) & !is.na(tdtrpnum)) %>%
  select(sampno,perno,tdtrpnum,whytrp1s) %>%
  rename(HOUSEID = sampno,
         PERSONID = perno,
         TDTRPNUM = tdtrpnum) %>%
  inner_join(triploc)
#923,572 trips
##### 3. match o-d tract to micro/geotypes #########


#tic()
#od.tr.purpose = trips %>%
#  mutate(o_geotype = match.geotype(o_geoid),
#         o_microtype = match.microtype(o_geoid),
#         d_geotype = match.geotype(d_geoid),
#         d_microtype = match.microtype(d_geoid))
#toc()

tic()
od.tr.purpose = trips %>%
  left_join(linker, by = c("o_geoid" = 'GEOID')) %>%
  rename(o_geotype = geotype,
         o_microtype = microtype) %>%
  left_join(linker, by = c("d_geoid" = "GEOID")) %>%
  rename(d_geotype = geotype,
         d_microtype = microtype)
toc()


tmp = as.data.frame(od.tr.purpose)

nhts.lookup.var<- function(var, table = NULL, year = 2017){
  # function to look up values of a variable in nhts data
  # note that 2009 data has no excel codebook, so there is not searchable dictionary to use
  if(year == 2017){
    load('./Data/nhts.codebook.17.RData')
    lookuptable = lookuptable17
  }else{
    stop('no year found')
  }
  if(var %in% unique(lookuptable$name)){
    if(is.null(table)){
      return(lookuptable[lookuptable$name == var,])
    }else{
      if(table %in% unique(lookuptable$table)){
        return(lookuptable[lookuptable$name == var & lookuptable$table == table,])
      }else{
        stop('no table found')
      }
    }
  }else{
    stop('no variable found')
  }
}
var = 'whytrp1s'
tmp.label = nhts.lookup.var(var)
if(class(tmp[,var]) %in% c('integer','numeric')){
  tmp[,var] = factor(tmp[,var],levels = as.numeric(tmp.label$value), labels = tmp.label$label)
}else{
  tmp[,var] = factor(tmp[,var],levels = tmp.label$value, labels = tmp.label$label)
}

od.tr.purpose = tmp


save(od.tr.purpose,file= file.path(datdir,'NHTS_tract_od_with_label.RData'))

od.tr.purpose2 = od.tr.purpose %>%
  select(-o_geoid,-d_geoid)

save(od.tr.purpose2,file= file.path(datdir,'NHTS_tract_od_with_label_no_location.RData'))


