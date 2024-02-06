* Merging all data for the FHWA geotypes
	* US CENSUS LODES
	* ACS 5 YR 2013-2017
	* CNT HTA Index
	* Slope data from USGS
	* Open Street Maps CBG Network Data
	* manipulated IPUMS Landcover summaries
	* commutes from LODES
	* climate data from NHGIS
	

* Set working directory
cd "C:/FHWA/For FHWA folks/microtype_input_preparation" 

global clean "CleanData"
global data "RawData"
global network "RawData"
global landcover "RawData/Landcover" 

* set output date
global outputdate "20191125"

* VARIABLES THAT NEED TO BE INCLUDED
* ROW
	* total right of way (by functional class?)
* network statistics
* population
* employment
* land area
* commute distances
* manufacturing employment 
* agricultural employment
* some new econ variables

**************
* LODES DATA FOR CROSSWALK AND EMPLOYMENT DENSITIES
************************
* Lodes data uses 2017 census tracts
import delimited $data/us_xwalk_tract_2017.csv, delimiter(comma) clear
*drop v1
rename cty fips_county
format tract %11.0f
rename tract trct

tempfile risso
save "`risso'"

****************
* merging with LEHD WAC data
***************
import delimited $data/wac_tract_2017.csv, delimiter(comma) clear
format trct %11.0f

merge 1:1 trct using `risso' 
* 506 obs don't match
drop _merge

tempfile harbor
save "`harbor'"

* NOTE: change this xwalk to 2015 bgrp level. 
* bring in LODES xwalk from 2015 to see if more match
import delimited $data/us_xwalk_cbg_2015.csv, delimiter(comma) clear
merge m:1 trct using `harbor'

* only 25 job tracts that still have no match for geography
drop _merge 
duplicates drop

* census tracts that end in 990000 to 990099 are tracts with  water areas
tostring trct , g(trct_str) format(%17.0g)

foreach x of varlist jobs* {
replace `x' = 0 if mi(`x') & strpos(trct_str, "990000")
}
drop trct_str jobs_ag
rename trct tract

* jobs
foreach x of varlist jobs_manu jobs_mi {
g pct_`x' = `x'/jobs_total 
drop `x'
}

label var pct_jobs_mi "Pct of jobs in NAICS 21 (tract)"
label var pct_jobs_manuf "Pct of jobs in NAICS 31-33 (tract)"
label var jobs_total "Total number of jobs (tract)" 
label var cbsa "Core-Based Statistical Area"
label var st_code "State acronym"
label var tract "Census tract"

tempfile right 
save "`right'" 

***********
* Adding slope data
********
import delimited $data/Slopes/slope_data.csv, delimiter(comma) clear varnames(1)

format bgrp %20.0f
tostring bgrp, g(bgrp_str) format (%20.0g)

bysort bgrp_str: egen slope_bgrp = mean(slope_mean)
drop bgrp_str slope_mean
duplicates drop

*duplicates drop

merge 1:m bgrp using `right'
drop _merge v1

bysort tract: egen slope_tract = mean(slope_bgrp)
drop slope_bgrp bgrp 
duplicates drop

tempfile gray
save "`gray'"

****************************
* THE FOLLOWING VARIABLES ARE ALL AT THE TRACT LEVEL
***************************

********************
* Importing ACS population data 
************
import delimited $data/acs_data_tracts.csv, delimiter(comma) clear		
keep geoid population 

*split geoid, p("US")
*drop geoid geoid1
*order geoid2
*destring(geoid2), replace
format geoid %11.0f
rename geoid tract 

merge 1:m tract using `gray'
* drop any tracts that are not from 2017
drop if _merge == 2
drop _merge 

tempfile blue
save "`blue'"

************
* Importing OSMnx network data
* Details and units: https://osmnx.readthedocs.io/en/stable/osmnx.html#module-osmnx.stats
* VARIABLES NEEDED: circuity_avg intersect_density_km street_length_total self_loop_proportion 
	*	degree_centrality_avg node_density_km street_length_avg dead_end_proportion
**************

import delimited "$data/usa-tracts-street_network-stats.csv", delimiter(comma) clear
rename state st_code
format geoid %20.0f
rename geoid tract 

keep area_km tract st_code circuity_avg street_density intersect_density_km street_length_total self_loop ///
		 street_length_avg dead_end_proportion

rename intersect_density_km intersection_density_km

rename area_km aland_km2017
label var aland_km2017 "land area of census tract in sq km"
label var street_length_avg "Avg street segment length (meters)"
label var circuity_avg "total edge length divided by sum of great circle distances between the nodes indecent to each edge"
label var self_loop_proportion "proportion of edges that have a single incident node"
label var intersection_density_km "intersection count divided by area (# per sq km)"
label var dead_end_proportion "proportion of nodes that are dead-ends"

g street_total_km = street_length_total /1000 
label var street_total_km "Total length of streets (km per tract)"

replace circuity = "." if circuity == "inf"
destring(circuity_avg), replace

* format decimal points 
foreach x of varlist street_length_avg {
format `x' %10.1f
}

merge 1:m tract using `blue' 
drop _merge 

* 422 tracts not matching. 

order tract st_code cbsa cbsaname ctyname 

tempfile bottlenose
save "`bottlenose'"

**************
* LAND USE SUMMARIES BY CENSUS TRACT
* Source: https://www.nhgis.org/user-resources/environmental-summaries
/*
	AREA:		total census tract area
	AREA11:    open water 
	AREA12:    perennial ice/snow
	AREA21:    developed, open space = less than 20% impervious surface
	AREA22:    developed, low intensity = 20-49% impervious surface
	AREA23:    developed, medium intensity = 50-79% impervious surfac
	AREA24:    developed, high intensity = 80% or more impervious surface
	AREA31:    barren land (rock/sand/clay)
	AREA41:    deciduous forest 
	AREA42:    evergreen forest 
	AREA43:    mixed forest 
	AREA52:    shrub/scrub
	AREA71:    grassland/herbaceous   
	AREA81:    pasture/hay
	AREA82:    cultivated crops 
	AREA90:    woody wetlands 
	AREA95:    emergent herbaceous wetlands */
************************************************************

import delimited $landcover/imp_surf_ag_tracts.csv, delimiter(comma) clear varnames(1)
drop lsad tractce name _imp_surfm statefp countyfp aff
g pct_water = awater/(awater + aland)
rename _ag_mean pct_ag_land
g dev_intensity = _imp/ 100
drop _imp

rename awater awater2017
rename aland aland2017

format geoid %13.0g
rename geoid tract

merge 1:m tract using `bottlenose'
drop _merge

tempfile northern 
save "`northern'"

* import land use for AK

import delimited $landcover/AK_ag_dev.csv, delimiter(comma) clear varnames(1)
g dev_intensity_ak = (_dev_open*.10) + (_low_dev*.35) + (_med_de * .65) + (_hi_dev * .90) 
g pct_water_ak = awater/(aland + awater)

rename aland aland2010
rename geoid tract

keep tract aland pct dev _ag_

merge 1:1 tract using `northern'

replace pct_ag_land = _ag_mean if _merge == 3 & mi(pct_ag_land)
replace dev_intensity = dev_intensity_ak if _merge == 3 & mi(dev_intensity)
replace aland2017 = aland2010 if _merge == 3 & mi(aland2017)
replace pct_water = pct_water_ak if _merge == 3 & mi(pct_water)
drop _merge _ag_mean dev_intensity_ak aland2010 pct_water_ak


tempfile wedell
save "`wedell'"

* import land use for Hawaii
import delimited $landcover/HI_ag_dev.csv, delimiter(comma) clear varnames(1)

foreach x of varlist _agmean _ag_mau _oahu_ag _mol_ag_me _nii_ag _lan_agm _hi_agmean {
replace _ag_mean = `x' if mi(_ag_mean)
drop `x'
}
foreach x of varlist _devmean_1 _maui_dev _oahu_de _mol_devme _nii_dev _hi_devmea _lan_dev {
replace _devmean = `x' if mi(_devmean)
drop `x'
}

rename geoid tract
rename aland aland2010
g pct_water_hi = awater/(aland + awater)

g dev_int_hi = _devmean*.9
drop _devmean awater

merge 1:1 tract using `wedell'

* NOTE: Hawaii only has open space and high development intensity. 
replace pct_ag_land = _ag_mean if _merge == 3 & mi(pct_ag_land)
replace dev_intensity = dev_int_hi if _merge == 3 & mi(dev_intensity)
replace aland2017 = aland2010 if _merge == 3 & mi(aland2017)
replace pct_water = pct_water_hi if _merge == 3 & mi(pct_water)

drop _merge dev_int_hi _ag_mean aland2010 pct_water_hi

*export dataset of census tract areas for commute calcs
preserve
keep tract aland2 awater
export  delimited $data/tract_areas_2017.csv, replace 
save $data/tract_areas_2017, replace
restore
*****************************
* IMPORTING HPMS and ARNOLD ROW BY CENSUS TRACT
******************************
merge m:1 tract using $data/row.dta
drop _merge 

forvalues i = 1/7 {
g pct_class`i' = lm_tract_fsys`i'/lm_all_tract
}

* for now drop actual mileage, for the clustering
drop lm_tract_f* aland awater
g area2017 = aland2017 + awater2017
label var area2017 "2017 Total area of census tract (sq meters)"
label var aland2017 "2017 Land area of census tract (sq meters)"
label var awater2017 "2017 Water area of census tract (sq meters)"

* Calculating lane-mile density over just the portion of tract that is land
g lane_miles_sqkm = lm_all_tract/aland_km

***********************
* CLEANING
********************
drop if mi(tract)

********
* CALCULATING DENSITIES (UPDATE THESE)
****************
order aland2017 aland_km
replace aland_km = aland2017/1000000 if mi(aland_km)
replace aland_km = aland2017/1000000 if mi(aland_km)

*drop aland2* 
g pop_density = population / aland_km
g job_density = jobs_total/ aland_km
label var pop_density "Population per sq km (tract)"
label var job_density "Jobs per sq km (tract)"

***********
* THE FOLLOWING VARIABLES ARE AT THE COUNTY OR CBSA LEVEL
*********************
rename cty county_fips
********************************
* Hi -TECH jobs from BLS --> Xiaodan's notes: cannot find the data
*******************************
*merge m:1 county_fips using $data/jobs_bls.dta
*drop if _merge == 2 
*drop _merge area_title*

tempfile stellar
save "`stellar'"

*****************
* ADDING AIR QUALITY DATA
******************
import delimited $data/air_quality_tracts.csv, delimiter(comma) clear
*drop v1
g pollutant_count = pb1978 + ozone2008 + so2 + ozone2015 + pm1997 + pm2006 + ///
					pm2012 + co + no2 + pm10 + so2_71 + pb2008
keep geoid poll
rename geoid tract
duplicates drop

merge 1:1 tract using `stellar' 
drop _merge 

tempfile beluga
save "`beluga'"

************
* ADDING ROAD GRADE DATA
***********
import delimited $data/road_grades_tract.csv, delimiter(comma) clear
merge 1:1 tract using `beluga'
drop _merge 


*replace grade = slope if mi(grade)

* using FHWA road grade if available
* if not, use road grade generated from graphopper commutes 
* if not, use slope data generated from natl elevation dataset
g grade_best = grade_mean
replace grade_best = grade if mi(grade_mean)
replace grade_best = slope if mi(grade_best)
count if mi(grade_best)

drop slope grade_mean grade
rename grade_best road_grade

tempfile narwal 
save "`narwal'"

*********
* ADDING BROADBAND CONNECTIONS
* Source link: https://www.fcc.gov/sites/default/files/tract_map_jun_2017.zip  
**************
import delimited $data/FCC/tract_map_jun_2017/tract_map_jun_2017.csv, delimiter(comma) clear
drop pcat_10
rename tractc tract
rename pcat_all broadband

merge 1:1 tract using `narwal'


drop _merge 

***********
* CLEANING
************
gen fips_state=state
drop if mi(county_fips) 
drop fips_st county_fips  state_code
replace fips_state = "15" if st_code == "HI"
replace fips_state = "02" if st_code == "AK"
replace pop_density = 0 if population == 0
g lanemeters_per_capita = lm_all_tract*1609.34/populat

* census tracts that end in 990000 to 990099 are tracts with  water areas
tostring tract , g(tract_str) format(%17.0g)
g water = 0
replace water = 1 if strpos(tract_str, "990000")
replace water = 1 if pct_water == 1

foreach x of varlist pct_jobs_m* {
replace `x' = 0 if job_density == 0
}

g jobs_housing_bal = jobs_total/population

g pct_hiway = pct_class1 + pct_class2
g pct_local_roads = pct_class6+ pct_class7
g pct_mid_roads = pct_class3+ pct_class4+ pct_class5
drop pct_class1-pct_class7

* for those missing HPMS data, replace with OSMNX data
replace lm_all_tract = street_total_km/1.609 if lm_all_tract == 0
replace lane_miles_sqkm = lm_all_tract/aland_km2017 if lane_miles_sqkm == 0

* assume two lanes
g street_meters_per_capita = 2*(street_total_km*1000)/population
replace lanemeters_p= street_meters_per_ca if lanemeters_per == 0
drop street_meter

* drop county-level variables from the first stage clusters
*drop jobs_hi_tech jobs_pct_

* almost half of census tracts are missing info for HPMS sample variables...
drop terrain signal pct_sigs_c 

preserve
* checking out how many non-cbsas are 
keep cbsa* ctyname fips* st_
duplicates drop
sort fips_c cbsa
bysort cbsa: gen dup = cond(_N==1, 0, _n)
restore

* drop un-needed variables for clustering
*drop jobs_total population lm* street_total state cbsaname ctyname area* aland* awater* tract_str street_length_total avg_speed 

drop jobs_total population lm* street_total state cbsaname ctyname area* aland* awater* tract_str street_length_total 

* number of missing variables by observation
egen hmiss = rowmiss(*)
tab hmiss
sort hmiss

*drop if hmiss >23 & water == 0

drop hmiss 

order tract st_code fips_state fips_county cbsa 
export delimited $clean/microtypes_inputs.csv, replace 






