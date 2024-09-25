#Creating user classes
#Variables
# Output file: 
# NATALIE POPOVICH
# NOV 11 2021
#############################

mywd <- "C:/FHWA/For FHWA folks/trip_generation_and_gems_inputs"
setwd(mywd)

# figuredir <- "Figures"
inputsdir <- "RawData"
datadir <- "CleanData"

library(tidyverse)
library(data.table)

##########
# LOAD DATA
###########3

# tract level income data from the NHTS
# NOTE: where is this file generated from?
inc <- fread(file.path(datadir, "NHTS/nhts_no_ids_1hrtimebins_with_imputation.csv")) %>%
      select(med_inc, HOUSEID) %>%
      unique()

#  nhts_no_ids_1hrtimebins_with_imputation.csv
# this file only has 217K obs 

# Raw (anonymized) data from NHTS
people <- fread(file.path(inputsdir, "NHTS-Public/perpub.csv")) %>%
  left_join(inc) # was merging by person ID and household. only need HOUSEID
#this file has 264K observations. we lose a bunch of folks in the merge
  

##############
# HHFAMINC   Household income bin
###############
#-9=Not ascertained
#-8=I don't know
#-7=I prefer not to answer
#01=Less than $10,000
#02=$10,000 to $14,999
#03=$15,000 to $24,999
#04=$25,000 to $34,999
#05=$35,000 to $49,999
#06=$50,000 to $74,999
#07=$75,000 to $99,999
#08=$100,000 to $124,999
#09=$125,000 to $149,999
#10=$150,000 to $199,999
#11=$200,000 or more

#################
# cleaning data for use in mode choice model
# Nov 11 2021 : Population groups include age, income, and vehicle ownership
##################
keep <- people %>% 
  select(HOUSEID, PERSONID, WTPERFIN,  
        HHFAMINC,  R_AGE_IMP, R_SEX_IMP, med_inc, HHVEHCNT) %>%
  unite('person_indx', HOUSEID:PERSONID, sep = "_", remove = F) %>%
  rename(age = R_AGE_IMP, income = HHFAMINC) %>%
  distinct() %>% 
  mutate(income = ifelse(income <1, NA, income))

#income 
sum(is.na(keep$income)) # 7975 obs

# impute income with tract-level median income if missing 
# assign tract-level median income to same scale as NHTS
keep <- keep %>% 
  mutate(tract_inc = case_when(
                        med_inc < 10000 ~ 1,
                        med_inc >= 10000 & med_inc < 15000 ~ 2, 
                        med_inc >= 15000 & med_inc < 25000 ~ 3, 
                        med_inc >= 25000 & med_inc < 35000 ~ 4, 
                        med_inc >= 35000 & med_inc < 50000 ~ 5, 
                        med_inc >= 50000 & med_inc < 75000 ~ 6, 
                        med_inc >= 75000 & med_inc < 100000 ~ 7, 
                        med_inc >= 100000 & med_inc < 125000 ~ 8, 
                        med_inc >= 125000 & med_inc < 150000 ~ 9, 
                        med_inc >= 150000 & med_inc < 200000 ~ 10,
                        med_inc >= 200000 ~ 11 )) %>%
  mutate(tract_inc = as.integer(tract_inc),
         income = as.integer(income))

keep = keep %>%
  dplyr::mutate(income = case_when(
    !is.na(income) ~ income, 
    TRUE ~ tract_inc))

sum(is.na(keep$income)) #1183 left 

# for those observations left without an income, use the median (6)
keep$income[is.na(keep$income)] <- 6

#binary indicator for income 
keep = keep %>% 
    mutate(income_over_med = case_when(
                        income > 5 ~ 1, 
                        income < 6 ~ 0),
                       income_below_med = case_when(
                         income < 6 ~ 1, 
                         income  > 5 ~ 0)) %>%
  select(-tract_inc, - med_inc) 

#binary indicator for age and sex
keep = keep %>% 
  mutate(under16 = case_when(
                age < 16 ~ 1,
                TRUE ~ 0 ),
        over16 = case_when(
                age >= 16 ~ 1, 
                TRUE  ~ 0 ), 
        male= case_when(
              R_SEX_IMP ==01  ~ 1,
              TRUE ~ 0)) %>%
  select(-R_SEX_IMP)


########################
# creating discrete user class bins
# groups of income and age
keep = keep %>%
  mutate(PopulationGroupID = case_when(
    income_over_med == 1 & HHVEHCNT > 0 & age >= 65 ~ "HighIncVehSenior",
    income_over_med == 0 & HHVEHCNT > 0 & age >= 65 ~ "LowIncVehSenior",
    income_over_med == 1 & HHVEHCNT == 0 & age >= 65 ~ "HighIncNoVehSenior",
    income_over_med == 0 & HHVEHCNT == 0 & age >= 65 ~ "LowIncNoVehSenior",
    income_over_med == 1 & HHVEHCNT == 0 & age < 65 ~ "HighIncNoVeh",
    income_over_med == 0 & HHVEHCNT == 0 & age < 65 ~ "LowIncNoVeh",
    income_over_med == 1 & HHVEHCNT > 0 & age < 65 ~ "HighIncVeh",
    income_over_med == 0 & HHVEHCNT > 0 & age < 65 ~ "LowIncVeh"))
  #select(-PopulationGroupID_3bin) # Xiaodan's notes: this variable is no longer available, guess Natalie dropped it at some point but didn't update it here

fwrite(keep, file.path(datadir, "nhts_user_classes_inc_veh_sr.csv"), row.names = F)












