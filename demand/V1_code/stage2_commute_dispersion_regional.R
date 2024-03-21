# CODE TO ESTIMATE THREE COMMUTE METRICS FOR THE SECOND STAGE CLUSTERS
  # METRIC 1: demand-supply burden at the microtype level
  # METRIC 2: normalized HHI index of distribution of PMT across tracts within each spatial unit
              #normalized by lane miles 
  # METRIC 3: bottleneck propensity at tract level, aggregated to microtype level distribution
  # METRIC 4: microtype distribution in each region 

library(data.table)
library(dplyr)
library(purrr)
library(tidyr)

mywd <- "C:/Users/viktoriia/Documents/microtype_input_preparation"
setwd(mywd)
###########################################################
# start with main census tracts, 73,056 
xwalk <- read.csv("./RawData/us_xwalk_tract_2017.csv", header = T)
xwalk <- xwalk[, -c(1)]
colnames(xwalk)[6] <- "GEOID"
keep <- c("GEOID", "cty", "cbsa")

#calculate number of each microtype per SPATIAL ID
type <- read.csv("./RawData/microtypes_output_clara6.csv", header = T)
all <- merge(type[,c('tract', 'cluster6')], xwalk[keep], by = "GEOID")

#creating unique identifier for second stage
all$spatial_id <- ifelse(all$cbsa <99999, all$cbsa, all$cty)

# import centricity measures
#NOTE 20200214: only 72,245 tracts are in this dataset
pmt <- read.csv("./RawData/centricity_all_v2.csv", header = T)
keep <- c("GEOID", "pmt_tract", "pmt_region")
all <- merge(all, pmt[keep], by = "GEOID", all.x = T)

#import lane miles per tract
lm <- read.csv("./RawData/row.csv", header = T)
lm <- lm[,c("tract", "lm_all_tract")]
colnames(lm)[1] <- "GEOID"
table(is.na(lm$lm_all_tract))

#merge trips with lane miles per tract
all <- merge(all, lm, by = "GEOID", all.x = T)
table(is.na(all$lm_all_tract))

# 504 tracts have missing lane-mile data. 
# for tracts that are all water, assign zero lane miles
all$lm_all_tract <- ifelse(is.na(all$lm_all_tract) & all$cluster6 == 7, 0 ,all$lm_all_tract) 

#for the others, impute using median for now (n = 266)
all$lm_all_tract <- ifelse(is.na(all$lm_all_tract), median(all$lm_all_tract, na.rm = T),all$lm_all_tract) 
table(is.na(all$lm_all_tract))

################
# METRIC 1: demand-supply burden at the microtype level
# variables needed: PMT per tract and lane-miles per tract. 
##################
#Steps:   1 - create unique ID for microtype-geotype pair
#         2 - aggregate tract-level PMT and LM to the microtype level in each spatial ID
#        3 - calculate ratio of total PMT to total LM per microtype
############################################
# IMPORT MICROTYPE OUTPUT:
#type <- read.csv("~/Box/FHWA/Data/CleanData/microtypes_output_clara6.csv", header = T)
#all <- merge(all, type[,c(7:8)], by = "GEOID")

# step 1: create unique ID for each microtype-geotype pair
all <- data.table(all)
all <- all[order(spatial_id, cluster6),]
all[, ID := .GRP, by = .(spatial_id, cluster6)]

# step 2: aggregate tract-level LM to the microtype level in each spatial ID
all <- all %>% group_by(ID) %>% mutate(lm_microtype = sum(lm_all_tract, na.rm = T))

# aggregate pmt to microtype level in each spatial ID
table(is.na(all$pmt_tract))

#811 tracts missing PMT,
# for tracts that are all water, assign zero person miles
all$pmt_tract <- ifelse(is.na(all$pmt_tract) & all$cluster6 == 7, 0 ,all$pmt_tract) 

# for the other 632 tracts... assign PMT = median for now... 
all$pmt_tract <- ifelse(is.na(all$pmt_tract), median(all$pmt_tract, na.rm = T),all$pmt_tract) 

all <- all %>% group_by(ID) %>% mutate(pmt_microtype = sum(pmt_tract, na.rm = T))
table(is.na(all$pmt_microtype))
# step 3: calculate ratio of total PMT to total LM per microtype
all$ds_burden <- all$pmt_microtype/all$lm_microtype
export <- all[,c("spatial_id", "cluster6", "ds_burden")]
export <- export %>% distinct()

# reshape wide by spatial id for export (NEEDS TO BE A DATAFRAME)
export <- reshape(as.data.frame(export),
                  timevar = "cluster6",
                  idvar = "spatial_id",
                  direction = "wide")

export[is.na(export)] <- 0

# export demand-supply burden at microtype level
#fwrite(export, file = "~/Box/FHWA/Data/CleanData/ds_burden.csv", row.names = F)


################
# METRIC 2: HHI index of distribution of PMT across tracts within each spatial unit
#normalized by lane miles 
# variables needed: PMT per tract and lane-miles per tract. 
# Steps:  1 - calculate tract level PMT per LM
#         2 - aggregate tract-level PMT/LM to spatial ID
#         3 - calculate tract level portion of PMT/LM to regional PMT/LM
#         4 - calcuate the HHI
#         5 - step 5: generate N = number of tracts in each region
#         6 - calcuate the normalized HHI
############################################
# step 1: PMT per lane mile in each tract
all$pmt_per_lm_tract <- all$pmt_tract/all$lm_all_tract
table(is.na(all$pmt_per_lm_tract))

# replace pmt_per_tract with missing if it is infinite (lane-miles of zero)
is.na(all) <- sapply(all, is.infinite)

#step 2 - aggregate tract-level PMT/LM to spatial ID and diivide
all <- all %>% group_by(spatial_id) %>% mutate(pmt_per_lm_region = sum(pmt_per_lm_tract, na.rm= T))

# step 3 - calculate tract level portion of PMT/LM to regional PMT/LM
all$s_i <- all$pmt_per_lm_tract/all$pmt_per_lm_region

# step 4: calculate HHI for each spatial ID
# H = sum_i (s_i squared) for each spatial unit
all$s_i_squared <- all$s_i * all$s_i

all <- all %>% group_by(spatial_id) %>% mutate(h = sum(s_i_squared, na.rm = T))
table(is.na(all$h))

# step 5: generate N = number of tracts in each region
all <- all %>% group_by(spatial_id) %>% add_tally()

# step 6 - calcuate the normalized HHI
# H* = (H-1)N / (1 - (1/N))
all$num <- all$h - (1/all$n)
all$den <- 1 - (1/all$n)
all$hhi_normalized <- all$num/all$den

# replace hhi_normalized = 1 if number of tracts = 1
all$hhi_normalized[which(all$n == 1)] <- 1

keep <- all[,c("spatial_id", "hhi_normalized")]
keep <- keep %>% distinct()

export <- merge(export, keep, by = "spatial_id") # here no HHIs are missing


################
# METRIC 3: bottleneck propensity at tract level, aggregated to microtype level distribution
# variables needed: PMT per tract and lane-miles per tract. 
##################
#Steps:   1 - start with pmt_per_lm_tract from metric (2)
#         2 - extract Xth percentile (start with 90) for each spatial ID
#         3 - generate pct of each microtype representation's from those tracts
############################################

# step 1 - start with pmt_per_lm_tract from metric (2)
# step 2 - extract Xth percentile (start with 90) for each spatial ID

#vector of percentiles to calculate
p <- c(0.75, .8, .85, 0.9)

#create a list of functions, with one for each quantile, using purrr::map and purrr::partial.
# assign names to each function (useful for the output of summarize) using purrr::set_names
p_names <- map_chr(p, ~paste0( "pctile", .x*100))
p_funs <- map(p, ~partial(quantile, probs = .x, na.rm = TRUE)) %>% 
  set_names(nm = p_names)

# calculate percentiles for tract level pmt/lm
ptiles <- all %>% 
  group_by(spatial_id) %>% 
  summarize_at(vars(pmt_per_lm_tract), funs(!!!p_funs)) # this function has been soft deprecated, but still works as of Feb 2020

# generate indicators for each tract if they are above the percentile
all <- merge(all, ptiles, by = "spatial_id")

all <- mutate(all, in90th = ifelse(pmt_per_lm_tract > pctile90, 1, 0)) 
all <- mutate(all, in85th = ifelse(pmt_per_lm_tract > pctile85, 1, 0))    
all <- mutate(all, in80th = ifelse(pmt_per_lm_tract > pctile80, 1, 0))        

# step  3 - generate pct of each microtype representation's from those tracts
# for each set of tracts in the percentile limit, calculate proportion of occurences of each microtype
# sum number of tracts of each microtype in percentile sample
# divide by total number of tracts in percentile sample

# Numerators: counts of occurence microtype for each sample 
n1 <- all %>% filter(in90th ==1) %>% group_by(spatial_id, cluster6) %>% count(in90th, name = "n")
n2 <- all %>% filter(in85th ==1) %>% group_by(spatial_id, cluster6) %>% count(in85th, name = "n")
n3 <- all %>% filter(in80th ==1) %>% group_by(spatial_id, cluster6) %>% count(in80th, name = "n")

# calculating proportional representation of each microtype 
dfs <- list(n1,n2,n3) #group data frames together

wide <- lapply(dfs, function(x) reshape(as.data.frame(x[,-c(3)]), #only keep necessary columns. reshape only works on data.frames
               timevar = "cluster6",
               idvar = "spatial_id",
               direction = "wide"))
rm(n1,n2,n3, dfs)
# replace all missing values with 0
wide <- lapply(wide, function(x) replace(x, is.na(x), 0))

# calculate total number of tracts in percentile sample
wide <- lapply(wide, function(x) {
            mutate(x, 
                  sample_n = rowSums(select(x, contains("n."))))
})

# calculate percent of representation 
#single df: n90 <- n90 %>% mutate_at(vars(n90_type.3:n90_type.7) , funs(P = ./n90$sumVar))
wide <- lapply(wide, function(x) {
  mutate_at(x,
            vars(n.3:n.7) , funs(pct_bottleneck = ./x$sample_n))
})

# keep only percentage variables
# single df: n90 <- n90 %>% mutate(sumVar = rowSums(select(., contains("type"))))
wide <- lapply(wide, function(x){
      select(x, matches("pct|spatial_id") )
})


# histogram of distribution for each cutoff
# for now use 85th%ile since it has more variation than 90th percentile
pct80 <- wide[[2]] # this dataset only has 2,044 observations... 

export <- merge(export, pct80, by = "spatial_id", all.x = T)

#####################
# MICROTYPE REPRESENTATION IN EACH REGION
###########################
#calculate number of each microtype per SPATIAL ID
#creating unique identifier for second stage
type$spatial_id <- ifelse(type$cbsa <99999, type$cbsa, type$fips_county)

# total number of tracts per cbsa
den <- type %>% group_by(spatial_id) %>% tally()

# count of tracts of each type per cbsa
num <- type %>% group_by(spatial_id, cluster6) %>% tally()

df <- merge(num,den, by = "spatial_id", all = T)
rm(den, count)

df$pct <- df$n.x/df$n.y

## Long to wide for number of tracts of each type
w <- reshape(df[,c(1,2,5)], 
             timevar = "cluster6",
             idvar = c("spatial_id"),
             direction = "wide")

# replace NAs with zeros 
w[is.na(w)] <- 0

export <- merge(export,w, by = "spatial_id", all = T)

fwrite(export, file = "./CleanData/regional_metrics_v2.csv", row.names = F)













