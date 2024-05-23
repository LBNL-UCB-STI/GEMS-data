*set the working directory
global rawdir "C:\FHWA\For FHWA folks\user_and_externality_cost\RawData"
global datadir "C:\FHWA\For FHWA folks\user_and_externality_cost\CleanData"
cd "$rawdir"

import excel "2019-Fare-Database.xlsx", sheet("fare") firstrow clear
rename City city
save "fare.dta", replace

import delimited "tract_place.csv", clear
joinby city using "fare.dta"
drop if fare==.
duplicates drop tractcode mode, force
tostring tractcode, replace format("%12.0f")
replace tractcode="0"+tractcode if statecode<10
gen county=substr(tractcode, 1, 5)
save "fare.dta", replace

*county
use "fare.dta", clear
collapse (mean) fare mini maxi, by(county mode)
save "county.dta", replace

*geo_micro
import delimited "tract_cluster_results_labeled.csv", clear
rename geoid tractcode
tostring tractcode, replace format("%12.0f")
replace tractcode="0"+tractcode if fips_st<10
keep tractcode geo_micro_type microtype geotype
merge 1:m tractcode using "fare.dta", keep(3) nogen
save "temp.dta", replace

collapse (mean) fare mini maxi, by(geo_micro_type mode)
save "geomicro.dta", replace

use "temp.dta", clear
collapse (mean) fare mini maxi, by(microtype mode)
save "micro.dta", replace

use "temp.dta", clear
collapse (mean) fare mini maxi, by(geotype mode)
save "geo.dta", replace

*merge them together
import delimited "tract_cluster_results_labeled.csv", clear
rename geoid tractcode
tostring tractcode, replace format("%12.0f")
expand 3
sort tractcode
by tractcode:  gen GroupID = _n
gen mode="bus" if GroupID==1
replace mode="rail_l" if GroupID==2
replace mode="rail_c" if GroupID==3
drop GroupID

merge 1:1 tractcode mode using "fare.dta"
save "temp.dta", replace

keep if _merge==3
save "merge1.dta", replace

use "temp.dta", clear
keep if _merge==1
replace county=substr(tractcode, 1, 5)
keep tractcode -state county
merge m:1 county mode using "county.dta"
save "temp.dta", replace

keep if _merge==3
save "merge2.dta", replace

use "temp.dta", clear
keep if _merge==1
keep tractcode -state county
merge m:1 geo_micro_type mode using "geomicro.dta"
save "temp.dta", replace

keep if _merge==3
save "merge3.dta", replace

use "temp.dta", clear
keep if _merge==1
keep tractcode -state county
merge m:1 microtype mode using "micro.dta"
save "temp.dta", replace

keep if _merge==3
save "merge4.dta", replace

use "temp.dta", clear
keep if _merge==1
keep tractcode -state county
merge m:1 geotype mode using "geo.dta"
keep if _merge==3
save "merge5.dta", replace


use "merge1.dta", clear
append using "merge2.dta"
append using "merge3.dta"
append using "merge4.dta"
append using "merge5.dta"

keep tractcode county state fare mini maxi
rename county countycode
export excel using "$datadir\transitfare_tract.xlsx", replace firstrow(variables)
 
erase "fare.dta"
erase "county.dta"
erase "geomicro.dta"
erase "geo.dta"
erase "micro.dta"
erase "temp.dta"
erase "merge1.dta"
erase "merge2.dta"
erase "merge3.dta"
erase "merge4.dta"
erase "merge5.dta"
 
 
 