# code to aggregate transit costs by mode and region

# NATALIE POPOVICH
# LAWRENCE BERKELEY NATIONAL LAB
# LAST UPDATED: APR 23 2020
#######################

# set working directory and sub-directories
mywd <- "C:/FHWA/For FHWA folks/Road_network_and_transit_cost_generation"
setwd(mywd)

#figuredir <- "./Figures"
rawdir <- "./RawData"
datadir <- "./CleanData"

# Install/load packages
install.packages('openxlsx')
library(openxlsx)
install.packages('tidyverse')
library(tidyverse)

######
# LOAD DATA
#######

# zipcode to county xwalk
# Source link: https://www.huduser.gov/portal/datasets/usps/ZIP_COUNTY_032020.xlsx
cty <- openxlsx::read.xlsx(file.path(rawdir,"ZIP_COUNTY_032020.xlsx"), sheet = 1) %>%
  select(ZIP, COUNTY)

# regional tract to CBSA crosswalk
xwalk <- read.csv(file.path(datadir, "us_xwalk_tract_2017_withID.csv")) %>%
  select(cbsa, cbsaname, cty, ctyname, spatial_id) %>% 
  distinct() %>%
  mutate(cty = sprintf("%05d",cty))

########################
# IMPORT CAPITAL EXPENSES BY AGENCY
####################
# All of these datasets can be searched for here: https://www.transit.dot.gov/ntd/ntd-data 
# Source link: https://www.transit.dot.gov/ntd/data-product/2018-capital-expenses 
# note these are only for existing services 
cap <- openxlsx::read.xlsx(file.path(rawdir,"./NTD/Capital_Expenses_2018.xlsx"), sheet = 3)

#remove spaces in variable names #Xiaodan's update: replace the old 'funs' statement as it is deprecated
df <- cap %>% 
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
ag <- openxlsx::read.xlsx(file.path(rawdir,"./NTD/2018 Agency Info.xlsx"), sheet = 1) %>% 
  rename_all(list(~str_replace(., " ", "."))) %>% 
  rename_all(list(~str_replace(., " ", "."))) %>% 
  rename_all(list(~str_replace(., " ", "."))) %>% 
  select(NTD.ID, Reporter.Type, City, State, Zip.Code, Service.Area.Sq.Miles, Population, UZA.Name, Sq.Miles) %>%
  distinct()

#################
# OPERATING EXPENSES BY MODE
# Source link: https://www.transit.dot.gov/ntd/data-product/2018-operating-expenses
######################
op <-openxlsx::read.xlsx(file.path(rawdir,"./NTD/Operating Expenses_2018.xlsx"), sheet = 1) %>% 
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
serv <- openxlsx::read.xlsx(file.path(rawdir,"./NTD/Service_2018.xlsx"), sheet = 3) %>%
  select(-contains("Question"), - contains("...")) %>% 
  rename_all(~ str_replace(., " ", ".")) %>% 
  rename_all(~ str_replace(., " ", ".")) %>% 
  rename_all(~ str_replace(., " ", ".")) %>% 
  rename_all(~ str_replace(., " ", ".")) %>% 
  rename_all(~ str_replace(., " ", ".")) %>%
rename(avg_speed = "Average.Speed.(mi/hr)", avg_trip_length = "Average.Passenger.Trip.Length.(mi)", 
         max_trains = "Max.Trains.in.Operation", pax_per_hr = "Passengers.per.Hour", TOS = "Type.of.Service",
         veh_revenue_hours = "Vehicle.Revenue.Hours", veh_revenue_miles = "Vehicle.Revenue.Miles" ) %>% 
  select(NTD.ID, Mode, max_trains, avg_speed, avg_trip_length, pax_per_hr, 
         veh_revenue_hours, veh_revenue_miles, Train.Miles, Train.Hours, 
               Unlinked.Passenger.Trips, Passenger.Miles, Directional.Route.Miles, TOS) %>%
  distinct()

###############
# TRACK AND ROW BY MODE
# Source link: https://www.transit.dot.gov/ntd/data-product/2018-track-and-roadway
###################
tr <- openxlsx::read.xlsx(file.path(rawdir,"./NTD/Track and Roadway_2018.xlsx"), sheet = 3) %>% # track 
  rename(track_miles_rev = Total.Revenue.Service, 
                    track_miles_nonrev = "Non-Revenue.Service", 
                    TOS = Type.Of.Service) %>% 
  select(NTD.ID, Mode, track_miles_rev, track_miles_nonrev, TOS)

tr <- tr %>% 
  mutate(NTD.ID = as.character(NTD.ID))

rd <- openxlsx::read.xlsx(file.path(rawdir,"./NTD/Track and Roadway_2018.xlsx"), sheet = 4) %>% # roadway 
rename(fixed_guideway_miles = Exclusive.Fixed.Guideway, 
                    busway_exclusive_miles = "Exclusive.High-Intensity.Busway",
                    control_access_miles = "Controlled.Access.High.Intensity.Busway.Or.HOV", 
                    row_total_miles =Total.Miles, 
                    TOS = Type.Of.Service) %>%
            select(NTD.ID, Mode, fixed_guideway_miles, busway_exclusive_miles, control_access_miles, row_total_miles, TOS)
rd <- rd %>% 
  mutate(NTD.ID = as.character(NTD.ID))
# Join all datasets together 
df = df %>%
  left_join(ag, by = "NTD.ID") %>% 
  filter(State != "PR" & State != "VI" & State !="GU" & State != "AS") %>%
  left_join(op, by = c("NTD.ID", "Mode", "TOS")) %>%
  left_join(serv, by = c("NTD.ID", "Mode", "TOS")) %>%
  left_join(tr, by = c("NTD.ID", "Mode", "TOS")) %>%
  left_join(rd, by = c("NTD.ID", "Mode", "TOS")) %>%
  mutate(ZIP = sprintf("%05d",Zip.Code)) %>% # add leading zeroes back to zipcode (fix to 5 characters)
  left_join(cty, by = "ZIP") %>% 
  rename(cty = COUNTY) %>%
  select(-NTD.ID, - Zip.Code, -Population, -Reporter.Type, - UZA.Name) %>%
  distinct() %>% 
  left_join(xwalk, by = "cty")

# there are lots of duplicates here
# drop observations where the same location is assigned to many different counties 
# but has the same VOMS, cap costs, VMT, mode, and 

df <- df %>% 
  distinct(Mode, Mode.VOMS, City, total_cap_exp, Service.Area.Sq.Miles, veh_revenue_hours, .keep_all = T)

#write.csv(df, file.path(datadir,"transit_costs_merged_062723.csv"), row.names = F)








