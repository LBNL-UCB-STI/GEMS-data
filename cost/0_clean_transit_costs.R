# code to aggregate transit costs by mode and region

#######################

# set working directory and sub-directories
gc()
mywd <- "C:/FHWA_R2/Cost"
setwd(mywd)

#figuredir <- "./Figures"
rawdir <- "./RawData"
datadir <- "./CleanData"

# Install/load packages
#install.packages('openxlsx')
library(openxlsx)
#install.packages('tidyverse')
library(tidyverse)
library(dplyr)
library(stringr)
######
# LOAD DATA
#######

# zipcode to county xwalk
# Source link: https://www.huduser.gov/portal/datasets/usps/ZIP_COUNTY_032020.xlsx
cty <- read.csv(file.path(datadir,"ZIP_COUNTY_LOOKUP_2023.csv")) %>%
  mutate(geoid = str_pad(as.character(geoid), 11, pad = "0"))%>% 
  mutate(county = str_sub(geoid, 1, 5), zip = sprintf("%05d", zip))%>% 
  select(zip, county) %>% # 188563 rows
  distinct() # 54300 rows

# regional tract to CBSA crosswalk
xwalk <- read.csv(file.path(datadir, "cleaned_lodes8_crosswalk_with_ID.csv")) %>%
  select(cbsa, cbsaname, cty, ctyname, spatial_id) %>% 
  mutate(cty = sprintf("%05d",cty)) %>%
distinct() 

########################
# IMPORT CAPITAL EXPENSES BY AGENCY
####################
# All of these datasets can be searched for here: https://www.transit.dot.gov/ntd/ntd-data 
# Source link: https://www.transit.dot.gov/ntd/data-product/2018-capital-expenses 
# note these are only for existing services 
capital_cost <- openxlsx::read.xlsx(file.path(rawdir,"./NTD/Capital_Expenses_2018.xlsx"), 
                           sheet = 'Total Capital Expenses by Mode')

#remove spaces in variable names 
#Xiaodan's update: replace the old 'funs' statement as it is deprecated
capital_cost <- capital_cost %>% 
  rename_all(~str_replace(., " ", ".")) %>% 
  rename_all(~str_replace(., " ", ".")) %>% 
  select(NTD.ID, Mode, TOS, Mode.VOMS, Total) %>% 
  rename(total_cap_exp = Total) %>% 
  # lots of duplicate entries with zero capital costs for the same agency
  filter(total_cap_exp > 0 ) %>% 
  distinct()

###############
# TRANSIT AGENCY INFO
# Source link: https://www.transit.dot.gov/ntd/data-product/2018-annual-database-agency-information
###############
agency <- openxlsx::read.xlsx(file.path(rawdir,"./NTD/2018 Agency Info.xlsx"), 
                          sheet = 'Agency Info') 
agency <- agency %>% 
  rename_all(list(~str_replace(., " ", "."))) %>% 
  rename_all(list(~str_replace(., " ", "."))) %>% 
  rename_all(list(~str_replace(., " ", "."))) %>% 
  select(NTD.ID, Reporter.Type, City, State, Zip.Code, Service.Area.Sq.Miles, Population, UZA.Name, Sq.Miles) %>%
  distinct()

#################
# OPERATING EXPENSES BY MODE
# Source link: https://www.transit.dot.gov/ntd/data-product/2018-operating-expenses
######################
op_cost <-openxlsx::read.xlsx(file.path(rawdir,"./NTD/Operating Expenses_2018.xlsx"), 
                              sheet = 'Operating Expenses') 

op_cost <- op_cost%>% 
  rename_all(list(~str_replace(., " ", "."))) %>% 
  rename_all(list(~str_replace(., " ", "."))) %>% 
  rename_all(list(~str_replace(., " ", "."))) %>% 
  filter(Operating.Expense.Type == "Total") %>% 
  rename(total_op_exp = Total.Operating.Expenses) %>%
  select(NTD.ID, Mode, total_op_exp, TOS) %>%
  distinct()

#################
# SERVICE ATTRIBUTES BY MODE
# Source link: https://www.transit.dot.gov/ntd/data-product/2018-service
######################
service <- openxlsx::read.xlsx(file.path(rawdir,"./NTD/Service_2018.xlsx"), 
                               sheet = 'Annual Service Data By Mode') 

service <- service %>%
  select(-contains("Question"), - contains("...")) %>% 
  rename_all(~ str_replace(., " ", ".")) %>% 
  rename_all(~ str_replace(., " ", ".")) %>% 
  rename_all(~ str_replace(., " ", ".")) %>% 
  rename_all(~ str_replace(., " ", ".")) %>% 
  rename_all(~ str_replace(., " ", ".")) %>%
rename(avg_speed = "Average.Speed.(mi/hr)", 
       avg_trip_length = "Average.Passenger.Trip.Length.(mi)",
       max_trains = "Max.Trains.in.Operation", 
       pax_per_hr = "Passengers.per.Hour", 
       TOS = "Type.of.Service",
       veh_revenue_hours = "Vehicle.Revenue.Hours", 
       veh_revenue_miles = "Vehicle.Revenue.Miles" ) %>% 
  select(NTD.ID, Mode, max_trains, avg_speed, avg_trip_length, pax_per_hr, 
         veh_revenue_hours, veh_revenue_miles, Train.Miles, Train.Hours,
         Unlinked.Passenger.Trips, Passenger.Miles, Directional.Route.Miles, TOS) %>%
  distinct()

###############
# TRACK AND ROW BY MODE
# Source link: https://www.transit.dot.gov/ntd/data-product/2018-track-and-roadway
###################
track <- openxlsx::read.xlsx(file.path(rawdir,"./NTD/Track and Roadway_2018.xlsx"), 
                             sheet = 'Track by Mode') 
track <- track %>% # track 
  rename(track_miles_rev = Total.Revenue.Service,
         track_miles_nonrev = "Non-Revenue.Service",
         TOS = Type.Of.Service) %>% 
  select(NTD.ID, Mode, track_miles_rev, track_miles_nonrev, TOS)

track <- track %>% 
  mutate(NTD.ID = as.character(NTD.ID))

roadway <- openxlsx::read.xlsx(file.path(rawdir,"./NTD/Track and Roadway_2018.xlsx"), 
                          sheet = 'Roadway by Mode')

roadway <- roadway %>% # roadway 
rename(fixed_guideway_miles = Exclusive.Fixed.Guideway,
       busway_exclusive_miles = "Exclusive.High-Intensity.Busway",
       control_access_miles = "Controlled.Access.High.Intensity.Busway.Or.HOV", 
       row_total_miles =Total.Miles,
       TOS = Type.Of.Service) %>%
  select(NTD.ID, Mode, fixed_guideway_miles, busway_exclusive_miles, control_access_miles, row_total_miles, TOS)

roadway <- roadway %>% 
  mutate(NTD.ID = as.character(NTD.ID))
# Join all datasets together

df <- capital_cost %>%
  left_join(agency, by = "NTD.ID") %>% 
  filter(State != "PR" & State != "VI" & State !="GU" & State != "AS") %>%
  left_join(op_cost, by = c("NTD.ID", "Mode", "TOS")) %>%
  left_join(service, by = c("NTD.ID", "Mode", "TOS")) %>%
  left_join(track, by = c("NTD.ID", "Mode", "TOS")) %>%
  left_join(roadway, by = c("NTD.ID", "Mode", "TOS")) %>%
  mutate(zip = sprintf("%05d",Zip.Code)) %>% # add leading zeroes back to zipcode (fix to 5 characters)
  left_join(cty, by = "zip", relationship = "many-to-many") %>% 
  rename(cty = county)

df <- df  %>%
  select(-Zip.Code, -Population, -Reporter.Type, - UZA.Name) %>%
  distinct() %>% 
  left_join(xwalk, by = "cty")
# 2685 rows

# there are lots of duplicates here
# drop observations where the same location is assigned to many different counties 
# but has the same VOMS, cap costs, VMT, mode, and 

df <- df %>% 
  distinct(NTD.ID, Mode, Mode.VOMS, City, total_cap_exp, Service.Area.Sq.Miles, veh_revenue_hours, .keep_all = T)
# 2017 rows
write.csv(df, file.path(datadir,"transit_system_cost.csv"), row.names = F)








