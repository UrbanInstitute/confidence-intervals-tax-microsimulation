global EXTRACT_DIR "K:\TPC\SCRATCH\For Lek\Sloan\data"

log using "$EXTRACT_DIR\bsw_gen.log", text replace
set linesize 255

version 14

display "$EXTRACT_DIR"

global num_reps = 201

* ==============================
* (1) Generate bootstrap weights
* ==============================

* calibrate_bsw calibrates each ECI_lvlgroup's total for every bootstrap to match the original data's total
capture program drop calibrate_bsw
program define calibrate_bsw
	local pweight_cur `1'

	total ones [pweight=`pweight_cur'], over(ECI_lvlgroup)

	foreach jj of numlist 0/11 {
		local jjp1 = `jj' + 1
		replace `pweight_cur' = `pweight_cur'*Totals_org[1,`jjp1']/_b[ones:`jj'] if (ECI_lvlgroup==`jj')
	}
	foreach jj of numlist 100/111 {
		local jjp1 = 12 + (`jj' - 100) + 1
		replace `pweight_cur' = `pweight_cur'*Totals_org[1,`jjp1']/_b[ones:`jj'] if (ECI_lvlgroup==`jj')
	}
end


use "$EXTRACT_DIR\bsw0.dta", clear

svyset [pweight=tabwt], strata(WSAMP)

generate double ones = 1
svy: total ones, over(ECI_lvlgroup)
matrix Totals_org = e(b)
matrix list Totals_org
*** display _b[ones:0]

bysort newseq: assert (_N == 1) /* Must sort the data to replicate the seed */
bsweights bsw, reps($num_reps) n(-1) seed(546841) double calibrate(calibrate_bsw @)

foreach jj of numlist 0/11 {
	summarize tabwt bsw* if (ECI_lvlgroup==`jj'), sep(0)
}
foreach jj of numlist 100/111 {
	summarize tabwt bsw* if (ECI_lvlgroup==`jj'), sep(0)
}

keep  newseq ECI_lvlgroup WSAMP bsw*
order newseq ECI_lvlgroup WSAMP bsw*
save "$EXTRACT_DIR\bsw${num_reps}.dta", replace

if ($num_reps >= 201) {
	* For running a sizable yet manageble number of bootstraps
	keep  newseq ECI_lvlgroup WSAMP bsw1-bsw201
	order newseq ECI_lvlgroup WSAMP bsw1-bsw201
	save "$EXTRACT_DIR\bsw201.dta", replace
}

if ($num_reps >= 21) {
	* For a test run (to check that distribution tables match the standard tables)
	keep  newseq ECI_lvlgroup WSAMP bsw1-bsw21
	order newseq ECI_lvlgroup WSAMP bsw1-bsw21
	save "$EXTRACT_DIR\bsw21.dta", replace
}

log close
