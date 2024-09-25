 * Code to clean the ATLAS database from FHWA
 * DATASETS
			* 2017 HPMS 
			* 2017 ARNOLD
			* 2015 HPMS SAMPLE 
* NATALIE POPOVICH
* BERKELEY LAB
* LAST UPDATED: SEP 30 2020 	
* Update Mar 23 2021: creating total network length variable per tract (not just lane miles)
* Updated Mar 26 2021: create total length variable by functional system  (not just lane miles)
**********************************

* Set working directory
cd "C:/FHWA/For FHWA folks/network_and_transition_matrix" 

global clean "CleanData"
global data "RawData"
global hpms "RawData/HPMS2017"
global arnold "RawData/ARNOLD2017"
global sample "RawData/hpms_sample_2015"
* set output date
global outputdate "20191115"
 

 **************************
 * CENSUS TRACT LEVEL ROW SUM STATS FROM STATE-LEVEL HPMS 2017 
 ********************
filelist, dir($hpms) pattern(*.txt)		 
tempfile files
save "`files'"
local obs = _N		

*read cleaned alabama data from the sample 
import delimited $hpms/al_hpms_sample_tract_intersect_2015_small.csv, clear
keep f_system through_la geoid aland awater length route_id psr ///
 state_code iri access_con aadt* speed*
 
*Manipulating raw data
forvalues i=1/`obs' {
    use "`files'" in `i', clear
    local f = dirname + "/" + filename
	display "`f'"
    import delimited "`f'" , delimiter(comma) clear varnames(1)
	keep f_system through_la geoid aland awater urban_code length route_id psr ///
 state_code iri access_con route_name aadt* speed*
	tempfile save`i'
    save "`save`i''", replace
	} 
	
*Appending files
use `save1', clear
forvalues i=2/`obs' {
    append using "`save`i''", force
}
		
tempfile hpms
save "`hpms'"


* NOTE: SOME ROUTE IDS ARE NUMERIC SOME STRING

****************
* CENSUS TRACT LEVEL ROW STATS FOR ARNOLD
************

*csv files
filelist, dir($arnold) pattern(*.csv)	 
tempfile files
save "`files'"
local obs = _N		

*Manipulating raw data
forvalues i=1/`obs' {
    use "`files'" in `i', clear
    local f = dirname + "/" + filename
	display "`f'"
    import delimited "`f'" , delimiter(comma) clear varnames(1)
	keep f_system through_la geoid aland awater urban_code length route_id psr ///
	state_code iri access_con route_name aadt* speed*
	tempfile save`i'
    save "`save`i''", replace
	} 
	
*Appending files
use `save1', clear
forvalues i=2/`obs' {
    append using "`save`i''", force
}

append using `hpms' , force

rename geoid tract
format tract %11.0f 
format aland %20.0f 
format awater %20.0f 
rename through_la lanes

order tract f_system length lanes

tab lanes f_sys

**********************
* data cleaning
* arc map naming is different than QGIS naming
****************************
replace aadt_singl = aadt_single_unit if mi(aadt_singl)
replace aadt_combi = aadt_combination if mi(aadt_combi)
replace lanes = through_lanes if mi(lanes)
replace access_con = access_control_ if mi(access_con)
replace speed_limit = speed_limi if mi(speed_limit)
drop aadt_single_unit aadt_combination through_lanes access_control_  speed_limi

* make assumptions for crazy number of lanes
* if lanes = 0 and f_sys  6 or 7, lanes = 2 for local roads (per Tianjia Tang)

replace lanes = 2 if mi(lanes) & inlist(f_sys, 6, 7)

* urban highways with the most lanes
* https://www.fhwa.dot.gov/policyinformation/tables/01.pdf
* Georgia ( I-75 - 15 lanes)
* California
replace lanes = 14 if lanes > 14 & urban_code == 51445 & strpos(route_id, "405") /* CA I-405 */
replace lanes = 12 if lanes > 14 & inlist(urban_code, 51445, 57709) & strpos(route_id, "005") /* CA I-5 */
replace lanes = 12 if lanes > 14 & inlist(urban_code,78904,90028) & strpos(route_id, "080") /* CA I-80 */
replace lanes = 12 if lanes > 14 & urban_code == 78661 & strpos(route_id, "015") /* CA I-15 */
* replace lanes = 8 if lanes == 0 & strpos(route_id, "210") & inlist(tract,6037101400,6037101300) /* CA I-210 */

* need to figure out naming convention for route_id. then could potentially clean data, for now leave as is.
* seems like everything with zero lanes ends with a _S suffix.


* Maryland
replace lanes = 14 if lanes > 14 & urban_code == 4843

* Wisconsin
replace lanes = 2 if lanes == 22 & tract == 55043961200

tab lanes f_sys

******************
* CLEANING UP ZERO-LANE OBSERVATIONS
*****************************

tab state if lanes == 0
* All arizona obs with zero lanes have an indicator in the route ID that it is different from the road segment
drop if state == 4 & lanes == 0 

* Virgina: combines East and Westbound routes and puts all lanes for one direction only
drop if state == 51 & lanes == 0 
drop if state == 19 & lanes == 0 /*Iowa does the same thing */

* MI: has duplicates for each road segment with iri = 0 and lanes = 0, drop all of them
drop if state == 27 & lanes == 0

* for now, drop all obs with zero lanes, seems like most of these are states that 
* separately measure each direction of traffic
drop if lanes == 0

drop urban_code

* need to aggregate for each system-lane # pair for each tract
* then reshape to wide so each tract has one observation 
g length_miles = length /1609.344
g lane_miles = length_miles*lanes

drop route_name length

***************
* INTERNATIONAL ROUGHNESS INDEX (IRI)
* available only for NHS roadways 
* and all segments of class 1, 2, and 3
**********************************

* aggregate all functional classes by IRI condition

* if IRI is missing, use psr
replace iri = . if iri == 99095
replace iri = . if iri == 0

******
* PSR
********
* where iri is missing, use PSR if available

* 4-5: new/almost new = 
* 3-4: few, if any, visible signs of surface deterioration.
* 2-3:  riding qualities of pavements in this category are noticeably 
* inferior to those of new pavements, and may be barely tolerable for high-speed traffic.
* 1-2: Pavements in this category have deteriorated to such an Extent that they 
* affect the speed of free-flow traffic. Flexible pavement may have large potholes and deep cracks.
* 0-1: Pavements in this category are in an extremely deteriorated condition. 
* The facility is passable only at reduced speeds, and with considerable ride discomfort.

* PSR 0-1: assign 99th% percentile IRI value = 377
replace iri = 377 if mi(iri) & psr < 1.0001
* PSR 1-2: assign 95th% percentile IRI = 250
replace iri = 250 if mi(iri) & psr >= 1 & psr < 2.0001
* PSR 2-3: assign "acceptable rating: 95 < IRI < 170
replace iri = 150 if mi(iri) & psr >= 2 & psr < 3.0001
* PSR 3-4: assign "good" rating 0 <IRI < 90
replace iri = 90 if mi(iri) & psr >= 3 & psr < 4.0001
* PSR 4-5: assign excellent rating percentile IRI = 10
replace iri = 10 if mi(iri) & psr >= 4 



* NUMERATOR: weighted IRI per lane mile
g iri_lm = iri*lane_miles
bysort tract: egen iri_num = sum(iri_lm)

* IRI DENOMINATOR: total lane miles of all classes per tract that have IRI value
bysort tract: egen lm_iri = sum(lane_miles) if !mi(iri)
bysort tract: egen iri_den = max(lm_iri)
g avg_iri = iri_num/iri_den

bysort tract: egen lm_all_tract = sum(lane_miles)
drop  iri* lm_iri psr

*aggregate by f_system-length
bysort tract f_sys: egen lm_tract_fsys = sum(lane_miles)

* total length of network for each functional system in miles
bysort tract f_sys: egen network_length_miles_fsys = sum(length_miles)
*  Xiaodan's notes variable rename from 'sys_length_miles_fsys' to 'network_length_miles_fsys' during bug fixing

* total network length in miles
* this should be the same as summing across all "network_length_miles_fsys" but it isn't...
bysort tract: egen sys_length_miles_tract = sum(length_miles) 

* this seems ways too high
*drop length_tract_sum route_id #Xiaodan notes: variable 'length_tract_sum' does not exist
drop route_id
************
* ACCESS CONTROL
*************

*Calculate percentage of total lane miles that are controlled access
bysort tract: g lm_full_control = sum(lane_miles) if access == 1
bysort tract: egen lm_full = max(lm_full_control)
replace lm_full = 0 if mi(lm_full)
bysort tract: g lm_part_control = sum(lane_miles) if access == 2
bysort tract: egen lm_part = max(lm_part)
replace lm_part = 0 if mi(lm_part)

g pct_controlf = lm_full/lm_all_tract
g pct_controlp = lm_part /lm_all_tract
drop lm_full* lm_part* access 

***************
* SPEED LIMITS: only required for some functional class systems
*****************
replace speed_limit = . if inlist(speed_limit, 0, 999)

* NUMERATOR: weighted speed limit per lane mile
g speed_lm = speed_limi*lane_miles
bysort tract: egen num = sum(speed_lm)

* DENOMINATOR: total lane miles of all classes per tract that have speed limit
bysort tract: egen den = sum(lane_miles) if !mi(speed_limi)
g avg_speed_limit = num/den

drop speed_lm num den speed_limit

foreach x of varlist avg_speed aadt_combi {
bysort tract: replace `x' = `x'[_n-1] if mi(`x')
bysort tract: replace `x' = `x'[_n+1] if mi(`x')
bysort tract: replace `x' = `x'[_n+1] if mi(`x')
bysort tract: replace `x' = `x'[_n+1] if mi(`x')
bysort tract: replace `x' = `x'[_n+1] if mi(`x')
bysort tract: replace `x' = `x'[_n+1] if mi(`x')
bysort tract: replace `x' = `x'[_n+1] if mi(`x')
bysort tract: replace `x' = `x'[_n+1] if mi(`x')
bysort tract: replace `x' = `x'[_n+1] if mi(`x')
bysort tract: replace `x' = `x'[_n+1] if mi(`x')
bysort tract: replace `x' = `x'[_n+1] if mi(`x')
bysort tract: replace `x' = `x'[_n+1] if mi(`x')
bysort tract: replace `x' = `x'[_n+1] if mi(`x')
}




************
* AADT
***********

*overall AADT (weighted by link length)
g aadt_link_length = aadt *length_miles
* I think this (bleow) is really VMT
bysort tract: egen aadt_tract = sum(aadt_link_length)



* CALCULATING FREIGHT TRAVEL IN TERMS OF VMT PER DAY (THEN PER HOUR)
*want: sum of freight AADT times link length in each census tract
g vmt_combi_link = aadt_combi *length_miles
bysort tract: egen vmt_combi_tract = sum(vmt_combi_link)

g vmt_sin_link = aadt_sin *length_miles
bysort tract: egen vmt_single_tract = sum(vmt_sin_link)



* percent of total AADT that is combination trucks
* DENOMINATOR:
gen aadt_lm = lane_miles * aadt
bysort tract: egen aadt_lm_tract = sum(aadt_lm)

*numerator: total combination AADT
gen combi_lm = lane_miles * aadt_combi
bysort tract: egen aadt_combi_tract = sum(combi_lm)
 
g pct_aadt_combi =  aadt_combi_tract /aadt_lm_tract

 * Truck AADT per lanemile 
*NUMERATOR: total combination AADT
* DENOMINATOR: total lane_miles in each tract
* NOTE AUG 13: I think this should actually just be by link length, not lane-link-length
bysort tract: egen den = sum(lane_miles) if !mi(aadt_combi)
g aadt_combi_per_lm = aadt_combi_tract /den

 
* total single truck AADT
*gen sin_lm = lane_miles * aadt_sin
*bysort tract: egen aadt_sin_tract = sum(sin_lm)
 
 **************
* reshape to wide so each tract has one observation 
**************
* need to drop obs that vary over census tracts in order to reshape
drop combi_lm den aadt_lm aadt_combi aadt_singl vmt_combi_link vmt_sin_link lane_miles aadt avg_speed length_miles lanes aadt_link_length

* unique id for each tract and functional class 
*egen id = group(tract f_sys)
duplicates drop
drop if mi(f_sys)
drop if f_sys == 0

* replace combi per lane mile with obsverations for the same tract
foreach x of varlist aadt_combi_per_lm {
bysort tract: replace `x'= `x'[_n+1] if mi(`x')
}

sort tract f_sys
quietly by tract f_sys: gen dup = cond(_N==1, 0, _n)
drop if dup > 0 & mi(aadt_combi_per_lm)
drop dup

bysort tract: replace aadt_combi_per_lm = aadt_combi_per_lm[_n+1] if mi(aadt_combi_per_lm)

duplicates drop

reshape wide lm_tract_fsys network_length_miles_fsys, i(tract) j(f_sys)

foreach x of varlist lm_tract_fsys* network_length_miles_fsys*{
replace `x'= 0 if mi(`x')
}

export delimited $clean/row_no_sample.csv, replace

tempfile row
save "`row'"
*************************************************************
**************
*IMPORTING HPMS 2015 FOR SAMPLE VARIABLES	
* road grade, terrain type, signal timing, 
********************
*************************************************
import delimited $sample/hpms_sample_lengths.csv, delimiter(comma) clear
tab state_code
*drop lsad affgeoid state_code tractce name countyfp lane_wid median_wid through_la aland awater #Xiaodan's note, keep last two variables for scaling
drop lsad affgeoid state_code tractce name countyfp lane_wid median_wid through_la

order sample_id begin_poin end_p
sort sample_id begin_poin
*drop length
duplicates drop

* if same segment has same starting and ending poitn but intersects with two different tracts, 
* drop second tract observation

sort sample_id begin_p end_p
quietly by sample_id begin_p end_p: gen dup = cond(_N==1, 0, _n)
* less than 2.5% of sample
drop if dup > 1 
drop dup

* create lenght for each sample ID, not portions of the sample
g segment_length_miles = end_poin - begin_poin

bysort sample_id: egen sample_length_miles = sum(segment_length)

rename geoid tract
format tract %11.0f 
format aland %20.0f 
format awater %20.0f 

order tract sample_id length segment_length sample_length state aland awater length	



********
* TERRAIN TYPE: ONLY REQUIRED FOR RURAL PLACES CLASSES 1-5
* * 1 - level, 2 - rolling, 3 - mountainous
**********
replace length = "0" if length == "nan"
destring length, replace
label var length "length of road segment in census tract i"

* terrain type, average by tracts
replace terrain_ty = . if terrain_ty == 0

bysort tract: egen terrain_den = sum(length) if !mi(terrain_ty)
g terrain_length = terrain_ty * length 
bysort tract: egen terrain_num = sum(terrain_length)
g terrain_intensity = terrain_num/terrain_den
drop terrain_length terrain_num terrain_den terrain_ty
replace terrain_int = 1 if terrain_in < 1


*********
* ROAD GRADE :* ONLY REQUIRED FOR CLASSES 1-4 (RURAL) AND 1-3 (URBAN)
* grades are in miles, length is in meters
* A 0 - .4
* B .5 - 2.4
* C 2.5 - 4.4 
* D 4.5 - 6.4
* E 6.5 - 8.4
* F 8.5 or greater
**********
* calculate percent of each segment of each grade
g ignore = 0

replace ignore = 1 if (grades_a == 0 & grades_b == 0 & grades_c == 0 & grades_d == 0 ///
			& grades_e == 0 & grades_f == 0)

g grades_length = grades_a + grades_b + grades_c + grades_d + grades_e + grades_f
replace grades_length = . if ignore == 1

* calculate percent of each grade of each road segment
foreach x of varlist grades* {
g pct_`x' = `x' /grades_length
}

* expand to each census tract length, the percentage of each grade
*denominator = total length of segment in each tract, only counted if have
* observable values grades
bysort tract: egen grades_den = sum(length) if ignore == 0 

* using median of each grade bin 
g grades_avg = .2*pct_grades_a + 1.45*pct_grades_b + 3.45*pct_grades_c + ///
					5.45*pct_grades_d + 7.45*pct_grades_e + 10*pct_grades_f
g grades_num = grades_avg * length if ignore == 0
bysort tract: egen grades_sum = sum(grades_num)
bysort tract: gen grade_mean = grades_sum/grades_de
drop pct_gra* ignore grades*
					
					
***********
* STOP SIGNS: * all sample roads except class rural 6 and 7 and urban 7
***********
* assumes that stop_sign density applies equally for all parts of the segment
g stop_density = stop_signs/segment_length

tab stop_signs
* stop signs look like there is complete data


*denominator = total length of segment in each tract, only counted if have
* observable values of stop_signs
bysort tract: egen den = sum(length) if !mi(stop_signs)

* numerator: weighted density  of stop signs by length
g stop_length = stop_density * length if !mi(stop_signs)
bysort tract: egen num = sum(stop_length)
bysort tract: gen stop_sign_density = num/den
drop num den stop_length stop_density stop_signs
label var stop_sign_density "Number of stop signs per mile (sample only)"


************
* SIGNAL TYPE:  all sample roads except class rural 6 and 7 and urban 7
* 1 Uncoordinated Fixed Time 
* 2 Uncoordinated Traffic Actuated.
* 3 Coordinated Progressive
* 4 Coordinated Real-time Adaptive
* 5 no signal
************
replace signal_ty = . if signal_ty == 0

* want proportion of signalized intersections that are coordinated
* denominator is the number of observations of road segments in each census tract
bysort tract: gen dup = cond(_N==1, 0, _n) if !mi(signal_ty)
bysort tract: egen den= max(dup)

* numerator: number of signalized intersections
g sig3 = 0
replace sig3 = 1 if signal_ty == 3
g sig4 = 0
replace sig4 = 1 if signal_ty == 4
bysort tract: egen sigs3t = sum(sig3)
bysort tract: egen sigs4t = sum(sig4)
g pct_sigs_coord = (sigs3t + sigs4t)/den
drop sig* dup den


***********
* NUMBER OF SIGNALIZED INTERSECTIONS: all sample roads except class rural 6 and 7 and urban 7
*************
* we don't know total number of intersections, so this is calculated per mile

* assumes that signal density applies equally for all parts of the segment
g sig_density = number_s/segment_length

*denominator = total length of segment in each tract, only counted if have
* observable values of signals
bysort tract: egen den = sum(length) if !mi(number_sig)

* numerator: weighted density  of stop signs by length
g sig_length = sig_density * length if !mi(number_sig)
bysort tract: egen num = sum(sig_length)
bysort tract: gen signal_density = num/den
drop num* den sig_*
label var signal_density "Number of signalized intersections per mile (sample only)"

sort sample_id tract

keep aland awater grade_mean pct_sigs terrain tract signal 
duplicates drop
sort tract

bysort tract: gen dup = cond(_N==1, 0, _n)
tab dup

foreach x of varlist terrain grade pct_sig signal {
bysort tract: replace `x' = `x'[_n-1] if mi(`x')
bysort tract: replace `x' = `x'[_n+1] if mi(`x')
bysort tract: replace `x' = `x'[_n+1] if mi(`x')
}
drop dup
duplicates drop


merge 1:1 tract using `row'


drop _merge	
	
save $clean/row.dta, replace

export delimited $clean/row.csv, replace

*export delimited $clean/row2.csv, replace


import delimited $clean/row.csv, clear

count if mi(lm_all_tract)


































