*set the working directory
global rawdir "C:\FHWA\For FHWA folks\user_and_externality_cost\RawData"
global datadir "C:\FHWA\For FHWA folks\user_and_externality_cost\CleanData"
cd "$rawdir"

import delimited "Uber_fare.csv", clear
save "Uber_fare.dta", replace

import delimited "tract_place.csv", clear
joinby city using "Uber_fare.dta"
tostring tractcode, replace format("%12.0f")
replace tractcode="0"+tractcode if statecode<10
gen county=substr(tractcode, 1, 5)
duplicates drop tractcode, force
save "Uber_fare.dta", replace

*county level average
use "Uber_fare.dta", clear
collapse (mean) minimumfare priceperminute pricepermile, by(county)
save "uber_county.dta", replace

*geo_micro level average
import delimited "tract_cluster_results_labeled.csv", clear
rename geoid tractcode
tostring tractcode, replace format("%12.0f")
replace tractcode="0"+tractcode if fips_st<10
keep tractcode geo_micro_type microtype geotype
merge 1:1 tractcode using "Uber_fare.dta", keep(3) nogen
save "uber_temp.dta", replace

use "uber_temp.dta", clear
collapse (mean) minimumfare priceperminute pricepermile, by(geo_micro_type)
save "uber_geomicro.dta", replace

use "uber_temp.dta", clear
collapse (mean) minimumfare priceperminute pricepermile, by(microtype)
save "uber_micro.dta", replace

use "uber_temp.dta", clear
collapse (mean) minimumfare priceperminute pricepermile, by(geotype)
save "uber_geo.dta", replace

*merge them together
import delimited "tract_cluster_results_labeled.csv", clear
rename geoid tractcode
tostring tractcode, replace format("%12.0f")
merge 1:1 tractcode using "Uber_fare.dta"
save "uber_temp.dta", replace

keep if _merge==3
save "uber_merge1.dta", replace

use "uber_temp.dta", clear
keep if _merge==1
replace county=substr(tractcode, 1, 5)
keep tractcode -state county
merge m:1 county using "uber_county.dta"
save "uber_temp.dta", replace

keep if _merge==3
save "uber_merge2.dta", replace

use "uber_temp.dta", clear
keep if _merge==1
keep tractcode -state county
merge m:1 geo_micro_type using "uber_geomicro.dta"
save "uber_temp.dta", replace

keep if _merge==3
save "uber_merge3.dta", replace

use "uber_temp.dta", clear
keep if _merge==1
keep tractcode -state county
merge m:1 microtype using "uber_micro.dta"
save "uber_merge4.dta", replace

use "uber_merge1.dta", clear
append using "uber_merge2.dta"
append using "uber_merge3.dta"
append using "uber_merge4.dta"

keep statecode tractcode county minimumfare priceperminute pricepermile
rename county countycode
export delimited using "$datadir\uber_tract_matched.csv", replace

*erase the temp datasets
erase "Uber_fare.dta"
erase "uber_county.dta"
erase "uber_geomicro.dta"
erase "uber_geo.dta"
erase "uber_micro.dta"
erase "uber_temp.dta"
erase "uber_merge1.dta"
erase "uber_merge2.dta"
erase "uber_merge3.dta"
erase "uber_merge4.dta"


