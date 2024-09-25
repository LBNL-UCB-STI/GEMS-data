* Created by Michelle Chen on 3/4/2023
* Modified by Hung-Chia Yang on 3/7/2023 to separate the R code from Stata do.file, and to update the directory path

*set the working directory
global dir "C:\FHWA\For FHWA folks\CleanData"
cd "$dir"
clear all

* auto
import delimited "parking_tract.csv", clear
rename tractcode GEOID
replace parking=parking/2
merge 1:1 GEOID using "micro_geo_results_with_imputation_volpe.dta", keep(3) nogen /*merge the geomicro_type*/
collapse (mean) parking, by(MicrotypeID)
gen Mode="auto"
rename parking PerEndCost
gen PerStartCost=0
gen PerMinuteCost=.
gen DailyCostPerVehicle=.
gen PerMileCost=.
save "auto_temp.dta", replace

* uber
import delimited "uber_tract_matched.csv", clear
rename tractcode GEOID
merge 1:1 GEOID using "micro_geo_results_with_imputation_volpe.dta", keep(3) nogen /*merge the geomicro_type*/
collapse (mean) minimumfare priceperminute pricepermile, by(MicrotypeID)
rename (minimumfare priceperminute pricepermile) (PerStartCost PerMinuteCost PerMileCost)
gen DailyCostPerVehicle=.
gen PerEndCost=.
gen Mode="ridehail"
save "ride_temp.dta", replace

* transit
import excel "transitfare_tract_final.xlsx", sheet("Sheet1") firstrow clear // this "transitfare_tract_final.xlsx" is idential to "transitfare_tract.xlsx" in Box
destring tractcode, replace
rename tractcode GEOID
merge m:1 GEOID using "micro_geo_results_with_imputation_volpe.dta", keep(3) nogen /*merge the geomicro_type*/
replace mode="rail" if mode=="rail_l" | mode=="rail_c"
collapse (mean) fare, by(MicrotypeID mode)
rename (mode fare)(Mode PerStartCost)
gen PerMinuteCost=.
gen DailyCostPerVehicle=.
gen PerMileCost=.
gen PerEndCost=.
save "transit_temp.dta", replace


clear all
use "auto_temp.dta", clear
append using "ride_temp.dta"
append using "transit_temp.dta"

export delimited using "modecost_03072023.csv", replace


erase "auto_temp.dta"
erase "ride_temp.dta"
erase "transit_temp.dta"

/*
* Move the following R script to a separate R code file  --> "C:\FHWA\For FHWA folks\Code\modecost_part2.R"
* =-=====================then use R to calculate fuel cost per mile: here are the r codes
# mywd <- "C:/Users/qmchen/Downloads"
setwd(mywd)
library(dplyr)


load("monetary_cost.RData")
load("travelmode.test.data.national_transgeo.RData")

geomicrotype <- TravelMode %>% distinct(trip_indx, o_microtype, o_geotype)
monetary_trip <- trips_veh %>% merge(geomicrotype, by="trip_indx")
modecost <- monetary_trip %>% group_by(o_microtype, o_geotype) %>%
  summarise(pricepermile=mean(2.41/mpg))
modecost <- modecost %>% filter(is.na(o_microtype)==F, is.na(o_geotype)==F)
modecost <- modecost %>% mutate(Microtype=paste0(o_geotype,"_", o_microtype))
*==============================================
* then manually copy the fuel cost to the csv permile for auto
*/