# Code to download Census LEHD data 
# aggregates raw data to census tract level for analysis
# We use only Workplace Area Characteristics and OD pairs
# https://lehd.ces.census.gov/data/
# exports tract-level data as .csv for clustering analysis
# Set working directory
mywd <- "C:/FHWA_R2/Demand"
setwd(mywd)
library(lehdr) # to download LEHD data from FTP
library(tidycensus) # to list state FIPS codes
library(data.table)
library(dplyr)



datadir <- "RawData"
cleandir <- "CleanData"
analysis_year <- 2021
###############
# DOWNLOAD LEHD DATA 
# NOTE: downloaded data are raw data. best to use imported data for analysis
#####################
us1 <- unique(fips_codes$state)[1:56]
# Workplace area characteristics (latest data from AK is 2016, and 2018 for AR and MS)
us <- us1[!us1 %in%  c("AS", "GU", "MP", "PR", "UM", "AK", "AR", "MS")]

# 2021 data for most us states
wac <- grab_lodes(us, 
                  analysis_year,
                  version = 'LODES8',
                  lodes_type = "wac",
                  job_type = "JT00", #all jobs combined
                  segment = "S000", # select total jobs
                  agg_geo = "tract")
head(wac)

# get 2016 data for AK
wac_ak <- grab_lodes('AK', 
                     2016,
                     version = 'LODES8',
                     lodes_type = "wac",
                     job_type = "JT00", #all jobs combined
                     segment = "S000", # select total jobs
                     agg_geo = "tract")

# get 2018 data for AK and MS
wac_ar_ms <- grab_lodes(c("AR", "MS"),
                     2018,
                     version = 'LODES8',
                     lodes_type = "wac",
                     job_type = "JT00", #all jobs combined
                     segment = "S000", # select total jobs
                     agg_geo = "tract")


#append AK, MS and AR to USA. LJ add: wac_us should be wac
wac_US <- rbind(wac, wac_ak, wac_ar_ms) %>% 
  select(year, state, w_tract, C000, CNS01, CNS02, CNS03, CNS04, CNS05, CNS06, CNS07, CNS08, CNS09, CNS10,
         CNS11, CNS12, CNS13, CNS14, CNS15, CNS16, CNS17, CNS18) 

colnames(wac_US) <- c('year', 'state','w_tract', 'total_jobs', "naics_11", "naics_21", "naics_22", "naics_23",
                      "naics_3133", "naics_42", "naics_4445", "naics_4849", "naics_51", "naics_52", "naics_53",
                      "naics_54", "naics_56", "naics_61", "naics_62", "naics_71", "naics_72")

#Export workplace area characteristics for each census tract
fwrite(wac_US, file = file.path(cleandir, paste0("wac_tract_", analysis_year, ".csv")), row.names = FALSE)

##################
# Origin-destination pairs (missing AK and SD for 2017)
#########################
for (st in us){
  print(st)
  ods.main <-grab_lodes(st, 
                        analysis_year,
                        version = 'LODES8',
                        lodes_type = "od",
                        job_type = "JT00", #all jobs combined
                        segment = "S000", # select total jobs
                        state_part = "main",
                        agg_geo = "tract") %>%
    select(year, state, w_tract, h_tract, S000)
  head(ods.main)
  
  ods.aux <- grab_lodes(st, 
                        analysis_year,
                        version = 'LODES8',
                        lodes_type = "od",
                        job_type = "JT00", #all jobs combined
                        segment = "S000", # select total jobs
                        state_part = "aux",
                        agg_geo = "tract")%>%
    select(year, state, w_tract, h_tract, S000) 
  
  ods <- rbind(ods.main, ods.aux)
  fwrite(ods, file.path(cleandir, "OD", paste0("od_pairs_", st, "_", analysis_year, ".csv")),  row.names = FALSE)
}



# Alaska from 2016

ods_ak.main  <- grab_lodes('AK', 
                           2016,
                           version = 'LODES8',
                           lodes_type = "od",
                           job_type = "JT00", #all jobs combined
                           segment = "S000", # select total jobs
                           state_part = "main",
                           agg_geo = "tract") %>%
  select(year, state, w_tract, h_tract, S000)

ods_ak.aux  <- grab_lodes('AK', 
                           2016,
                           version = 'LODES8',
                           lodes_type = "od",
                           job_type = "JT00", #all jobs combined
                           segment = "S000", # select total jobs
                           state_part = "aux",
                           agg_geo = "tract") %>%
  select(year, state, w_tract, h_tract, S000)

ods_ak <- rbind(ods_ak.main, ods_ak.aux)
fwrite(ods_ak, file.path(cleandir, "OD", "od_pairs_AK_2016.csv"),  row.names = FALSE)
# get 2018 data for AK and MS
ods_ar_ms.main  <- grab_lodes(c('AR', 'MS'), 
                                  2018,
                                  version = 'LODES8',
                                  lodes_type = "od",
                                  job_type = "JT00", #all jobs combined
                                  segment = "S000", # select total jobs
                                  state_part = "main",
                                  agg_geo = "tract") %>%
  select(year, state, w_tract, h_tract, S000)

ods_ar_ms.aux  <- grab_lodes(c('AR', 'MS'), 
                                2018,
                                version = 'LODES8',
                                lodes_type = "od",
                                job_type = "JT00", #all jobs combined
                                segment = "S000", # select total jobs
                                state_part = "aux",
                                agg_geo = "tract") %>%
  select(year, state, w_tract, h_tract, S000)

ods_ar_ms <- rbind(ods_ar_ms.main,ods_ar_ms.aux)
fwrite(ods_ar_ms, file.path(cleandir, "OD", "od_pairs_ARMS_2016.csv"),  row.names = FALSE)

############################
# CROSSWALK FOR COUNTY, STATE, AND CBSA  --> XXu note: this part has been performed under 0_clean_boundaries.R now
################################
#Note that the CBSAs change periodically and this crosswalk could be updated even without updating the LODES7 commute data

# xwalk <- fread(file.path(datadir, "LEHD/LODES2017/us_xwalk.csv")) %>%
#   select(st, stusps, stname, cty, ctyname, trct, cbsa, cbsaname) %>%
#   rename(fips_st = st, 
#        st_code = stusps, 
#        state = stname, 
#        tract = trct) %>% 
#   filter(fips_st < 60) %>% # remove all the US territories
#   distinct() # drop all duplicates since df was at block level
# 
# fwrite(xwalk, file.path(cleandir, "us_xwalk_tract_2017.csv"),  row.names = FALSE)
# #XIAODAN'S NOTES: this output may very likely be the same as us_xwalk_tract_2017_withID.CSV, probably just got slightly edited by Natalie
# 
# # 2015 crosswalk for census tracts that have changed
# xwalk <- fread(file.path(datadir, "LEHD/LODES2015/us_xwalk2015.csv")) %>%
#   select(st, stusps, stname, cty, ctyname, trct, cbsa, cbsaname) %>%
#   rename(fips_st = st,
#          st_code = stusps,
#          state = stname,
#          tract = trct) %>%
#   filter(fips_st < 60) %>% # remove all the US territories
#   distinct() # drop all duplicates since df was at block level
# fwrite(xwalk, file.path(cleandir, "us_xwalk_cbg_2015.csv"),  row.names = FALSE)


