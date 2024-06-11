*set the working directory
global rawdir "C:\FHWA\For FHWA folks\user_and_externality_cost\RawData"
global datadir "C:\FHWA\For FHWA folks\user_and_externality_cost\CleanData"
cd "$rawdir"

import excel "parkopedia.xlsx", clear sheet("Sheet1") firstrow
save "parking.dta", replace

import delimited "tract_place.csv", clear
joinby city using "parking.dta"
duplicates drop
keep tractcode parking
save "parking.dta", replace

import delimited "tract_place.csv", clear
merge 1:1 tractcode using "parking.dta"
replace parking=0 if parking==.
drop _merge

export delimited using "$datadir\parking_tract.csv", replace