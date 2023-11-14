*set the working directory
global rawdir "C:\FHWA\For FHWA folks\user_and_externality_cost\RawData"
global datadir "C:\FHWA\For FHWA folks\user_and_externality_cost\CleanData"
cd "$rawdir"

import delimited "survey_vehicles2019v", clear
duplicates drop make model, force
keep make model vehtype
rename make makecode
rename model modelcode
save "temp.dta", replace

import excel "vehmodel.xlsx", clear sheet("make") firstrow
save "make.dta", replace 
import excel "vehmodel.xlsx", clear sheet("model") firstrow
save "model.dta", replace 

use "temp.dta", clear
merge m:1 makecode using "make.dta", nogen
merge m:1 modelcode using "model.dta", nogen
duplicates drop make model, force
save "temp.dta", replace

import delimited "vehicles", clear
drop if year>2017
keep make model city08
joinby make model using "temp.dta"
rename city08 mpg 
save "temp0.dta", replace

collapse (mean) mpg, by(make model)
save "merge1.dta", replace

use "temp0.dta", clear
collapse (mean) mpg, by(make vehtype)
save "maketype.dta", replace

use "temp0.dta", clear
collapse (mean) mpg, by(vehtype)
save "type.dta", replace

use "temp0.dta", clear
collapse (mean) mpg, by(make)
save "make.dta", replace

*merge
use "temp.dta", clear
merge 1:1 make model using "merge1.dta"
save "temp0.dta", replace

keep if _merge==3
save "merge1.dta", replace

use "temp0.dta", clear
keep if _merge==1
drop mpg _merge
merge m:1 make vehtype using "maketype.dta"
save "temp0.dta", replace

keep if _merge==3
save "merge2.dta", replace

use "temp0.dta", clear
keep if _merge==1
drop mpg _merge
merge m:1 vehtype using "type.dta"
save "temp0.dta", replace

keep if _merge==3
save "merge3.dta", replace

use "temp0.dta", clear
keep if _merge==1
drop mpg _merge
merge m:1 make using "make.dta"
save "temp0.dta", replace

keep if _merge==3
save "merge4.dta", replace

use "temp0.dta", clear
keep if _merge==1
drop mpg _merge
*other left are motorcycle
*manually change the mpg
gen mpg=40
save "merge5.dta", replace


use "merge1.dta", clear
append using "merge2.dta"
append using "merge3.dta"
append using "merge4.dta"
append using "merge5.dta"
drop _merge
export delimited using "$datadir\veh_mpg.csv", replace

erase make.dta
erase maketype.dta
erase model.dta
erase type.dta
erase temp.dta
erase temp0.dta
erase merge1.dta
erase merge2.dta
erase merge3.dta
erase merge4.dta
erase merge5.dta

