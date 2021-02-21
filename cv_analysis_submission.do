// BP

/* load dataset: patient_dataextraction.dta (pre-analysis in Python) */
use patient_dataextraction.dta

/* MERGE DATASETS */
cd "~/Desktop"
merge 1:m patientunitstayid using bp_median, nogen // bp_median pre-analysis in Python
sort patientunitstayid hoursinceadm

/* POST MERGING PROCESSING */

	// drop if incomplete data
	drop if apachescore < 0 /* delete observation if apachescore <0 (which is not possible) */
	drop if apachescore == . | gender == . | age == . | mech_vent == . | sedatives == . | iculosdays == . | hoslosdays == . | mortality == . | vasodilators == . | vaso_ino == . | b_blockers == . | nonsel_ccb == . | ace_arb == . | sel_ccb == . | diuretics == .

	// drop if < 72 hours of data
	by patientunitstayid: egen max_hour = max(hoursinceadm)
	drop if max_hour < 72
	
	/* identify hour of admission for each patient */
	bys patientunitstayid (hoursinceadm): gen hourofadm = hourofday[1]

	/* count hourofday continuously */
	gen hourofdaycont = .
	replace hourofdaycont = hourofadm + hoursinceadm
	
	// observations / day
	bys patientunitstayid (hoursinceadm): egen obs_day1 = count(gender) if hoursinceadm <= 24
	by patientunitstayid (hoursinceadm): egen obs_day2 = count(gender) if hoursinceadm > 24 & hoursinceadm <= 48
	by patientunitstayid (hoursinceadm): egen obs_day3 = count(gender) if hoursinceadm > 48 & hoursinceadm <= 72
	by patientunitstayid (hoursinceadm): egen obs_day4 = count(gender) if hoursinceadm > 72 & hoursinceadm <= 96
	by patientunitstayid (hoursinceadm): egen obs_day5 = count(gender) if hoursinceadm > 96 & hoursinceadm <= 120
	by patientunitstayid (hoursinceadm): egen obs_day6 = count(gender) if hoursinceadm > 120 & hoursinceadm <= 144
	by patientunitstayid (hoursinceadm): egen obs_day7 = count(gender) if hoursinceadm > 144 & hoursinceadm <= 168
	by patientunitstayid (hoursinceadm): egen obs_day8 = count(gender) if hoursinceadm > 168 & hoursinceadm <= 192
	by patientunitstayid (hoursinceadm): egen obs_day9 = count(gender) if hoursinceadm > 192 & hoursinceadm <= 216
	by patientunitstayid (hoursinceadm): egen obs_day10 = count(gender) if hoursinceadm > 216 & hoursinceadm <= 240
	by patientunitstayid (hoursinceadm): egen obs_day11 = count(gender) if hoursinceadm > 240 & hoursinceadm <= 264
	by patientunitstayid (hoursinceadm): egen obs_day12 = count(gender) if hoursinceadm > 264 & hoursinceadm <= 288
	by patientunitstayid (hoursinceadm): egen obs_day13 = count(gender) if hoursinceadm > 288 & hoursinceadm <= 312
	by patientunitstayid (hoursinceadm): egen obs_day14 = count(gender) if hoursinceadm > 312 & hoursinceadm <= 336
	by patientunitstayid (hoursinceadm): egen obs_day15 = count(gender) if hoursinceadm > 336 & hoursinceadm <= 360
	
	by patientunitstayid (hoursinceadm): egen observations_day1 = max(obs_day1)
	by patientunitstayid (hoursinceadm): egen observations_day2 = max(obs_day2)
	by patientunitstayid (hoursinceadm): egen observations_day3 = max(obs_day3)
	by patientunitstayid (hoursinceadm): egen observations_day4 = max(obs_day4)
	by patientunitstayid (hoursinceadm): egen observations_day5 = max(obs_day5)
	by patientunitstayid (hoursinceadm): egen observations_day6 = max(obs_day6)
	by patientunitstayid (hoursinceadm): egen observations_day7 = max(obs_day7)
	by patientunitstayid (hoursinceadm): egen observations_day8 = max(obs_day8)
	by patientunitstayid (hoursinceadm): egen observations_day9 = max(obs_day9)
	by patientunitstayid (hoursinceadm): egen observations_day10 = max(obs_day10)
	by patientunitstayid (hoursinceadm): egen observations_day11 = max(obs_day11)
	by patientunitstayid (hoursinceadm): egen observations_day12 = max(obs_day12)
	by patientunitstayid (hoursinceadm): egen observations_day13 = max(obs_day13)
	by patientunitstayid (hoursinceadm): egen observations_day14 = max(obs_day14)
	by patientunitstayid (hoursinceadm): egen observations_day15 = max(obs_day15)
	
	// de-string admission time
	gen admtime = clock(unitadmittime24, "hms")
	format admtime %tc
	gen h = hh(admtime)
	summ h


// (linearized) single cosinor to test for significant 24 hour rhythm
gen sintime = sin(2*_pi*hourofdaycont/24)
gen costime = cos(2*_pi*hourofdaycont/24)

// quadratic term of hoursinceadm
gen hoursinceadmsq = hoursinceadm * hoursinceadm

// test signifiance of amplitude of every single cosinor

gen sig_single_cos_1 = .
gen sig_single_cos_2 = .
gen sig_single_cos_3 = .
gen sig_single_cos_4 = .
gen sig_single_cos_5 = .

gen amp_single_cos_1 = .
gen amp_single_cos_2 = .
gen amp_single_cos_3 = .
gen amp_single_cos_4 = .
gen amp_single_cos_5 = .

gen phase_single_cos_1 = .
gen phase_single_cos_2 = .
gen phase_single_cos_3 = .
gen phase_single_cos_4 = .
gen phase_single_cos_5 = .

// identify individual patients
egen x = group(patientunitstayid) 

cd "~/Desktop"
save "bp_data_pre_analyzed.dta", replace

clear



// single cosinor, BP

local ij 2500 5000 7500 10000 12500 15000 17500 20000 22500 25000

forvalues k = 1(1)4 {

	foreach var2 of local ij {

		cd "~/Desktop"
		use bp_data_pre_analyzed.dta 
		
		keep bp_median costime sintime x hoursinceadm observations_day2 observations_day3 observations_day4 observations_day5 observations_day6
		gen sig_single_cos_`k' = .
		gen amp_single_cos_`k' = .
		gen phase_single_cos_`k' = .
		
		tempvar ab1
		tempvar ab2
		tempvar ab3
		tempvar ab4

		gen `ab1' = 0
		gen `ab2' = 0
		gen `ab3' = 0
		gen `ab4' = 0

		replace `ab1' = 1 if hoursinceadm > 24 & hoursinceadm < 72 & observations_day2 > 5 & observations_day3 >5 & observations_day2 != . & 		observations_day3 != .

		replace `ab2' = 1 if hoursinceadm > 24 & hoursinceadm < 96 & observations_day2 > 5 & observations_day3 >5 & observations_day4 > 5 & 		observations_day2 != . & observations_day3 != . & observations_day4 != . 

		replace `ab3' = 1 if hoursinceadm > 24 & hoursinceadm < 120 & observations_day2 > 5 & observations_day3 >5 & observations_day4 > 5 & 		observations_day5 >5 & observations_day2 != . & observations_day3 != . & observations_day4 != . & observations_day5 != .

		replace `ab4' = 1 if hoursinceadm > 24 & hoursinceadm < 144 & observations_day2 > 5 & observations_day3 >5 & observations_day4 > 5 & 		observations_day5 >5 & observations_day6 > 5 & observations_day2 != . & observations_day3 != . & observations_day4 != . & 					observations_day5 != . & observations_day6 != . 
		
		keep if `ab`k'' == 1
	
		keep if x > (`var2' - 2500) & x <= `var2'
		
		summ x
		local min = r(min)
		local max = r(max)

		local ov bp_median
	
		forvalues i = `min'/`max' { 

			capture regress `ov' costime sintime if x == `i' 
			capture nlcom sqrt((_b[costime])^2+(_b[sintime])^2)
			capture matrix b = r(b)
			capture matrix V = r(V)	
			capture local b = (b[1,1])
			capture local se = sqrt(V[1,1])
			capture local ll = `b' - 1.96*`se'
			replace sig_single_cos_`k' = 1 if x == `i' & `ll' > 0 & `ll' != .
			capture noisily replace sig_single_cos_`k' = 0 if x == `i' & `ll' <= 0 
			capture noisily replace amp_single_cos_`k' = `b' if x == `i'

			capture if _b[sintime] >= 0 & _b[costime] >= 0 {
			capture nlcom (atan(_b[sintime]/_b[costime]) * (24/(2*_pi)))
			}
			capture else if _b[sintime] < 0 & _b[costime] >= 0 {
			capture nlcom ((atan(_b[sintime]/_b[costime]) * (24/(2*_pi))) + 24)
			}
			capture else if _b[sintime] >= 0 & _b[costime] < 0 {
			capture nlcom ((atan(_b[sintime]/_b[costime]) * (24/(2*_pi))) + 12)
			}
			capture else if _b[sintime] < 0 & _b[costime] < 0 {
			capture nlcom ((atan(_b[sintime]/_b[costime]) * (24/(2*_pi))) + 12)
			}
			capture matrix c = r(b)
			capture local c = (c[1,1])
			capture noisily replace phase_single_cos_`k' = `c' if x == `i'
		
			}
			
		cd "~/Desktop"
		save "cv_cos_`var2'_`k'_bp_v3.dta", replace
		clear

		}
		
	}
	
clear


// merge datasets

// single cosinor, BP

local ij 2500 5000 7500 10000 12500 15000 17500 20000 22500 25000

forvalues k = 1(1)4 {

	foreach var2 of local ij {
	
		cd "~/Desktop"
		use cv_cos_`var2'_`k'_bp_v3.dta
		
		bys x (hoursinceadm): keep if _n == 1
	
		keep x sig_single_cos_`k' amp_single_cos_`k' phase_single_cos_`k'
		
		local var sig_single_cos_`k' amp_single_cos_`k' phase_single_cos_`k'
		
		foreach i of local var {
			rename `i' `i'_bp
			}
	
		save "cv_cos_`var2'_`k'_bp_for_merge_v1.dta", replace
		clear
		
		}
	
	}
	

// append

local ik 5000 7500 10000 12500 15000 17500 20000 22500 25000

forvalues k = 1(1)4 {

	cd "~/Desktop"
	use cv_cos_2500_`k'_bp_for_merge_v1.dta

	foreach var3 of local ik {
	
		append using cv_cos_`var3'_`k'_bp_for_merge_v1.dta
		
		}
		
	save "cv_cos_all_bp_appended_`k'.dta", replace
	clear

	}
	
cd "~/Desktop"
use cv_cos_all_bp_appended_1.dta
merge 1:1 x using cv_cos_all_bp_appended_2.dta, nogen
merge 1:1 x using cv_cos_all_bp_appended_3.dta, nogen
merge 1:1 x using cv_cos_all_bp_appended_4.dta, nogen
save "cv_cos_all_bp_merged.dta", replace
clear




// combine BP and covariates

cd "~/Desktop"
use bp_data_pre_analyzed.dta 
bys x (hoursinceadm): keep if _n == 1
keep x patientunitstayid
save "x_pt_id_bp.dta", replace
clear

cd "~/Desktop"
use x_pt_id_bp.dta
merge 1:1 x using cv_cos_all_bp_merged.dta, nogen
drop x
save "cv_cos_bp_merged_pt_id.dta", replace
clear


cd "~/Desktop"
use patient_dataextraction.dta
cd "~/Desktop"
merge 1:1 patientunitstayid using cv_cos_bp_merged_pt_id.dta, nogen
keep if sig_single_cos_1_bp != . 

// combine non-sel and selective CCB
gen ccb = 0
replace ccb = 1 if nonsel_ccb == 1 | sel_ccb == 1

save "cos_for_analysis.dta", replace





// Analysis

// Association with baseline clinical characteristics

cd "~/Desktop"
putexcel set "table1_ind", modify

local v E
local w F
local x G
local y H

local b 1

local cov gender age apachescore explicit_sepsis organ_dysfunction mech_vent vasodilators vaso_ino sedatives b_blockers ccb ace_arb diuretics

foreach var of local cov {
	regress phase_23_mm_bp `var' 
	
	matrix results = r(table)
	putexcel `v'`b' = matrix(results[1,2])
	putexcel `w'`b' = matrix(results[5,2])
	putexcel `x'`b' = matrix(results[6,2])
	putexcel `y'`b' = matrix(results[4,1])
	local b = `b' + 1
	
	lincom _b[_cons] + _b[`var']
	local coef = r(estimate)
	local se = r(se)
	local p = r(p)
	putexcel `v'`b' = `coef'
	putexcel `w'`b' = (`coef' - (`se'*1.95))
	putexcel `x'`b' = (`coef' + (`se'*1.95))
	local b = `b' + 2
	}



	
	// age
	
// rescale age
summ age
replace age = age - (r(mean)-40)
replace age = age / 10

// regression

cd "~/Desktop"
putexcel set "table1_ind", modify

local v A
local w B
local x C
local y D

local b 1

regress amp_23_mm_bp age 
	matrix results = r(table)
	putexcel `v'`b' = matrix(results[1,2])
	putexcel `w'`b' = matrix(results[5,2])
	putexcel `x'`b' = matrix(results[6,2])
	putexcel `y'`b' = matrix(results[4,1])
	local b = `b' + 1
	
	lincom _b[_cons] + _b[age]
	local coef = r(estimate)
	local se = r(se)
	local p = r(p)
	putexcel `v'`b' = `coef'
	putexcel `w'`b' = (`coef' - (`se'*1.95))
	putexcel `x'`b' = (`coef' + (`se'*1.95))
	local b = `b' + 1
	
	lincom _b[_cons] + (2*_b[age])
	local coef = r(estimate)
	local se = r(se)
	local p = r(p)
	putexcel `v'`b' = `coef'
	putexcel `w'`b' = (`coef' - (`se'*1.95))
	putexcel `x'`b' = (`coef' + (`se'*1.95))
	local b = `b' + 1
	
	lincom _b[_cons] + (3*_b[age])
	local coef = r(estimate)
	local se = r(se)
	local p = r(p)
	putexcel `v'`b' = `coef'
	putexcel `w'`b' = (`coef' - (`se'*1.95))
	putexcel `x'`b' = (`coef' + (`se'*1.95))
	local b = `b' + 1
	
	lincom _b[_cons] + (4*_b[age])
	local coef = r(estimate)
	local se = r(se)
	local p = r(p)
	putexcel `v'`b' = `coef'
	putexcel `w'`b' = (`coef' - (`se'*1.95))
	putexcel `x'`b' = (`coef' + (`se'*1.95))
	local b = `b' + 1

	
	
	// apache score
	
// divide apachescore into quartiles
xtile apache_quartiles = apachescore, nq(4)

// regression

cd "~/Desktop"
putexcel set "table1_ind", modify

local v A
local w B
local x C
local y D

local b 1

regress amp_23_mm_bp i.apache_quartiles 
	matrix results = r(table)
	putexcel `v'`b' = matrix(results[1,5])
	putexcel `w'`b' = matrix(results[5,5])
	putexcel `x'`b' = matrix(results[6,5])
	local b = `b' + 1

	lincom _b[_cons] + (_b[2.apache_quartiles])
	local coef = r(estimate)
	local se = r(se)
	local p = r(p)
	putexcel `v'`b' = `coef'
	putexcel `w'`b' = (`coef' - (`se'*1.95))
	putexcel `x'`b' = (`coef' + (`se'*1.95))
	putexcel `y'`b' = matrix(results[4,2])
	local b = `b' + 1
	
	lincom _b[_cons] + (_b[3.apache_quartiles])
	local coef = r(estimate)
	local se = r(se)
	local p = r(p)
	putexcel `v'`b' = `coef'
	putexcel `w'`b' = (`coef' - (`se'*1.95))
	putexcel `x'`b' = (`coef' + (`se'*1.95))
	putexcel `y'`b' = matrix(results[4,3])
	local b = `b' + 1
	
	lincom _b[_cons] + (_b[4.apache_quartiles])
	local coef = r(estimate)
	local se = r(se)
	local p = r(p)
	putexcel `v'`b' = `coef'
	putexcel `w'`b' = (`coef' - (`se'*1.95))
	putexcel `x'`b' = (`coef' + (`se'*1.95))
	putexcel `y'`b' = matrix(results[4,4])
	local b = `b' + 1
	

// frequencies

cd "~/Desktop"
putexcel set "table1_freq", modify

local v A
local w B
local x C
local y D

local a 1
local b 2

local cov gender explicit_sepsis organ_dysfunction mech_vent vasodilators vaso_ino sedatives b_blockers ccb ace_arb diuretics

foreach var of local cov {

	tab `var' if amp_single_cos_1_bp != . & sig_single_cos_1_bp == 1, matcell(freq)
	
	putexcel `v'`a'=(freq[1,1])
	putexcel `w'`a'=(((freq[1,1])/r(N))*100)
	
	local a = `a' + 1
	
	putexcel `v'`a'=(freq[2,1])
	putexcel `w'`a'=(((freq[2,1])/r(N))*100)
	
	local a = `a' + 2
	
	}	
	
	

// testing significance of predictors

	// amplitude

local var amp_single_cos_1_bp

local a 1
local b 2

local cov apachescore age gender explicit_sepsis organ_dysfunction mech_vent vasodilators vaso_ino sedatives b_blockers ccb ace_arb diuretics

regress iculosdays c.`var'##c.`var' `cov' if mortality == 0 
estimates store m`a'
regress iculosdays `cov' if mortality == 0 & `var' != . 
estimates store m`b'
lrtest m`a' m`b'

local a 1
local b 2

local cov apachescore age gender explicit_sepsis organ_dysfunction mech_vent vasodilators vaso_ino sedatives b_blockers ccb ace_arb diuretics

regress hoslosdays c.`var'##c.`var' `cov' if mortality == 0 
estimates store m`a'
regress hoslosdays `cov' if mortality == 0 & `var' != . 
estimates store m`b'
lrtest m`a' m`b'

local a 1
local b 2

local cov apachescore age gender explicit_sepsis organ_dysfunction mech_vent vasodilators vaso_ino sedatives b_blockers ccb ace_arb diuretics

logistic mortality c.`var'##c.`var' `cov'  
estimates store m`a'
logistic mortality `cov' if `var' != . 
estimates store m`b'
lrtest m`a' m`b'

	// phase
	
local k bp

foreach var of local k {

	forvalues i=1(1)4 {
		gen phase_cat_`i'_`var' = .
		replace phase_cat_`i'_`var' = 0 if phase_single_cos_`i'_`var' >= 0 & phase_single_cos_`i'_`var' < 4
		replace phase_cat_`i'_`var' = 1 if phase_single_cos_`i'_`var' >= 4 & phase_single_cos_`i'_`var' < 8 
		replace phase_cat_`i'_`var' = 2 if phase_single_cos_`i'_`var' >= 8 & phase_single_cos_`i'_`var' < 12
		replace phase_cat_`i'_`var' = 3 if phase_single_cos_`i'_`var' >= 12 & phase_single_cos_`i'_`var' < 16 
		replace phase_cat_`i'_`var' = 4 if phase_single_cos_`i'_`var' >= 16 & phase_single_cos_`i'_`var' < 20 
		replace phase_cat_`i'_`var' = 5 if phase_single_cos_`i'_`var' >= 20 & phase_single_cos_`i'_`var' < 24 
		}
	
	}
	
	

local var phase_cat_4_bp

local a 1
local b 2

local cov apachescore age gender explicit_sepsis organ_dysfunction mech_vent vasodilators vaso_ino sedatives b_blockers ccb ace_arb diuretics

regress iculosdays i.`var' `cov' if mortality == 0 
estimates store m`a'
regress iculosdays `cov' if mortality == 0 & `var' != . 
estimates store m`b'
lrtest m`a' m`b'

local a 1
local b 2

local cov apachescore age gender explicit_sepsis organ_dysfunction mech_vent vasodilators vaso_ino sedatives b_blockers ccb ace_arb diuretics

regress hoslosdays i.`var' `cov' if mortality == 0
estimates store m`a'
regress hoslosdays `cov' if mortality == 0 & `var' != . 
estimates store m`b'
lrtest m`a' m`b'

local a 1
local b 2

local cov apachescore age gender explicit_sepsis organ_dysfunction mech_vent vasodilators vaso_ino sedatives b_blockers ccb ace_arb diuretics

logistic mortality i.`var' `cov'  
estimates store m`a'
logistic mortality `cov' if `var' != . 
estimates store m`b'
lrtest m`a' m`b'


// Outcome plots

	// Amplitude

regress iculosdays c.amp_single_cos_1_bp##c.amp_single_cos_1_bp apachescore age gender explicit_sepsis organ_dysfunction mech_vent vasodilators vaso_ino sedatives b_blockers ccb ace_arb diuretics if mortality == 0 

margins, at(amp_single_cos_1_bp=(3(0.1)15))
marginsplot, xlabel(3(0.1)15) recast(line) recastci(rarea)











	

