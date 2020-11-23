* [bs_dist_prog.do] Stata program for producing bootstraped Tax Model distributional estimates of the proposal $PROP_CUR

log using "$EXTRACT_DIR\\bs_dist_${PROP_CUR}${num_reps}.log", text replace

version 14

* The values of globals EXTRACT_DIR (data directory), num_reps (number of bootstrap repetitions) and PROP_CUR (current proposal) are defined in "main_reps##.do"
display "$EXTRACT_DIR"
display $num_reps
display "$PROP_CUR"

* +++++++++++++++++++++++++++++++++++++
* (1) Commands for capturing statistics
* +++++++++++++++++++++++++++++++++++++

* summarize_statorg_capture: Capture statistics without running bootstraps after running summarize
capture program drop summarize_statorg_capture
program define summarize_statorg_capture
	local matrow_num `1' /* Matrix row number */
	local matcol_num `2' /* Matrix column number */

	matrix storedorg_numobs[`matrow_num',`matcol_num'] = r(N)
	matrix storedorg_mean[`matrow_num',`matcol_num'] = r(mean)
	matrix storedorg_variance[`matrow_num',`matcol_num'] = r(Var)
	matrix storedorg_wgtobs[`matrow_num',`matcol_num'] = r(sum_w)
end

* ratio_statorg_capture: Capture statistics without running bootstraps after running ratio
capture program drop ratio_statorg_capture
program define ratio_statorg_capture
	local matrow_num `1' /* Matrix row number */
	local matcol_num `2' /* Matrix column number */

	matrix storedorg_numobs[`matrow_num',`matcol_num'] = e(_N)
	matrix storedorg_mean[`matrow_num',`matcol_num'] = e(b)
	matrix storedorg_variance[`matrow_num',`matcol_num'] = e(V)
	matrix storedorg_wgtobs[`matrow_num',`matcol_num'] = e(N)
end

* statbs_capture: Capture statistics after running bootstraps
capture program drop statbs_capture
program define statbs_capture
	local matrow_num `1' /* Matrix row number */
	local matcol_num `2' /* Matrix column number */

	matrix storedbs_numobs[`matrow_num',`matcol_num'] = e(N)
	matrix storedbs_mean[`matrow_num',`matcol_num'] = e(b)
	matrix storedbs_variance[`matrow_num',`matcol_num'] = e(V)
	matrix storedbs_bsmean[`matrow_num',`matcol_num'] = e(b_bs)
	matrix storedbs_bias[`matrow_num',`matcol_num'] = e(bias)
	matrix storedbs_z0[`matrow_num',`matcol_num'] = e(z0)
	matrix storedbs_bsse[`matrow_num',`matcol_num'] = e(se)
		matrix tmp_ci = e(ci_normal)
	matrix storedbs_lbcinormal[`matrow_num',`matcol_num'] = tmp_ci[1,1]
	matrix storedbs_ubcinormal[`matrow_num',`matcol_num'] = tmp_ci[2,1]
		matrix tmp_ci = e(ci_percentile)
	matrix storedbs_lbcip[`matrow_num',`matcol_num'] = tmp_ci[1,1]
	matrix storedbs_ubcip[`matrow_num',`matcol_num'] = tmp_ci[2,1]
		matrix tmp_ci = e(ci_bc)
	matrix storedbs_lbcibc[`matrow_num',`matcol_num'] = tmp_ci[1,1]
	matrix storedbs_ubcibc[`matrow_num',`matcol_num'] = tmp_ci[2,1]
end

* +++++++++++++++++++++++++++++++++++++++++++
* (2) Prepare the data for running bootstraps
* +++++++++++++++++++++++++++++++++++++++++++

use "$EXTRACT_DIR\baseline_2019cl.dta", clear
keep newseq depind tabwt eci_cl cl_burden
rename eci_cl cl_eci_cl

* Import the proposal's variables
merge newseq using "$EXTRACT_DIR\\${PROP_CUR}.dta", keep(alt_burden)
	assert (_merge == 3)
drop _merge
bysort newseq: assert (_N == 1)

* Import bootstrap weights
merge newseq using "$EXTRACT_DIR\bsw${num_reps}.dta"
	assert (_merge == 3)
drop _merge
bysort newseq: assert (_N == 1)

generate double ones = 1

generate double cl_AfterTaxIncome = cl_eci_cl - cl_burden

generate double FedTaxChange = alt_burden - cl_burden
foreach JJ of numlist 0/11 {
	generate double FedTaxChange`JJ' = FedTaxChange*(ECI_lvlgroup == `JJ')
}

generate double Winner = (FedTaxChange <= -10)
generate double Winner_amt = FedTaxChange*Winner
generate double Loser  = (FedTaxChange >=  10)
generate double Loser_amt  = FedTaxChange*Loser

generate double diff_AfterTaxIncome = -1*FedTaxChange

keep if (depind == 0) /* For distributional estimates, drop dependent filers */
drop depind
bysort newseq: assert (_N == 1)

* ECI-level group numbers to be used below 
matrix group_num = J(13,1,99)
foreach JJ of numlist 0/11 {
	local JJp1 = `JJ' + 1
	matrix group_num[`JJp1',1] = `JJ'
}
matrix list group_num

* +++++++++++++++++++++++++++++++++++++++++++++++++++++++
* (3) Replicate the proposal's summary distribution table
* +++++++++++++++++++++++++++++++++++++++++++++++++++++++

* Matrix for storing bootstrap estimates
* 12 ECI groups + 1 total = 13 rows, 9 variables

matrix storedorg_numobs   = J(13,1,-3), group_num, J(13,9,0) /* number of observations */
matrix storedorg_mean     = J(13,1,-2), group_num, J(13,9,0) /* avearge, original */
matrix storedorg_variance = J(13,1,-1), group_num, J(13,9,0) /* variance */
matrix storedorg_wgtobs   = J(13,1,-4), group_num, J(13,9,0) /* number of observations */

* --------------------
* Each ECI level group
* --------------------

foreach JJ of numlist 0/11 {
	local JJp1 = `JJ' + 1

	count if (ECI_lvlgroup == `JJ') & (Winner == 1)
	if (r(N) > 0) {
		* (1A) Number of Winners
		ratio Winner/ones [iweight=tabwt] if (ECI_lvlgroup == `JJ')
		ratio_statorg_capture `JJp1' 3

		* (1B) Avearge Tax Change among Winners
		*Q* Is this the right statistics for bootstraping since # winners vary across simulations???
		summarize Winner_amt [iweight=tabwt] if (ECI_lvlgroup == `JJ') & (Winner == 1)
		summarize_statorg_capture `JJp1' 4
	}
	else {
		summarize tabwt if (ECI_lvlgroup == `JJ')

		* (1A) Number of Winners
		local col_num = 3
		matrix storedorg_numobs[`JJp1',`col_num'] = r(N)
		matrix storedorg_mean[`JJp1',`col_num'] = 0
		matrix storedorg_variance[`JJp1',`col_num'] = 0
		matrix storedorg_wgtobs[`JJp1',`col_num'] = r(sum)


		* (1B) Avearge Tax Change among Winners
		local col_num = 4
		matrix storedorg_numobs[`JJp1',`col_num'] = r(N)
		matrix storedorg_mean[`JJp1',`col_num'] = 0
		matrix storedorg_variance[`JJp1',`col_num'] = 0
		matrix storedorg_wgtobs[`JJp1',`col_num'] = r(sum)
	}

	count if (ECI_lvlgroup == `JJ') & (Loser == 1)
	if (r(N) > 0) {
		* (1C) Number of Losers
		ratio Loser/ones [iweight=tabwt] if (ECI_lvlgroup == `JJ')
		ratio_statorg_capture `JJp1' 5

		* (1D) Avearge Tax Change among Losers
		*Q* Is this the right statistics for bootstraping since # losers vary across simulations???
		summarize Loser_amt [iweight=tabwt] if (ECI_lvlgroup == `JJ') & (Loser == 1)
		summarize_statorg_capture `JJp1' 6
	}
	else {
		summarize tabwt if (ECI_lvlgroup == `JJ')

		* (1C) Number of Losers
		local col_num = 5
		matrix storedorg_numobs[`JJp1',`col_num'] = r(N)
		matrix storedorg_mean[`JJp1',`col_num'] = 0
		matrix storedorg_variance[`JJp1',`col_num'] = 0
		matrix storedorg_wgtobs[`JJp1',`col_num'] = r(sum)

		* (1D) Avearge Tax Change among Loser
		local col_num = 6
		matrix storedorg_numobs[`JJp1',`col_num'] = r(N)
		matrix storedorg_mean[`JJp1',`col_num'] = 0
		matrix storedorg_variance[`JJp1',`col_num'] = 0
		matrix storedorg_wgtobs[`JJp1',`col_num'] = r(sum)
	}

	* (2) Percent Change in After-Tax Income = Total of diff_AfterTaxIncome / Total of cl_AfterTaxIncome
	ratio diff_AfterTaxIncome/cl_AfterTaxIncome [iweight=tabwt] if (ECI_lvlgroup == `JJ')
	ratio_statorg_capture `JJp1' 7

	* (3) Share of Total Federal Tax Change of group JJ = Total FedTaxChange`JJ' / Total FedTaxChange
	*T* For a revenue neutral proposal, FedTaxChange is approximately $0 and thus this ratio is not well-defined!
	ratio FedTaxChange`JJ'/FedTaxChange [iweight=tabwt]
	ratio_statorg_capture `JJp1' 8

	* (4) Average Federal Tax Change
	summarize FedTaxChange [iweight=tabwt] if (ECI_lvlgroup == `JJ')
	summarize_statorg_capture `JJp1' 9

	* (5) Average Federal Tax Rate Change (% points) = Total FedTaxChange / Total cl_eci_cl
	ratio FedTaxChange/cl_eci_cl [iweight=tabwt] if (ECI_lvlgroup == `JJ')
	ratio_statorg_capture `JJp1' 10

	* (6) Average Federal Tax Rate Under the Proposal = Total alt_burden / Total cl_eci_cl
	ratio alt_burden/cl_eci_cl [iweight=tabwt] if (ECI_lvlgroup == `JJ')
	ratio_statorg_capture `JJp1' 11
}

* -------------
* All tax units
* -------------

	local JJp1 = `JJp1' + 1

	count if (Winner == 1)
	if (r(N) > 0) {
		* (1A) Number of Winners
		ratio Winner/ones [iweight=tabwt]
		ratio_statorg_capture `JJp1' 3

		* (1B) Avearge Tax Change among Winners
		*Q* Is this the right statistics for bootstraping since # winners vary across simulations???
		summarize Winner_amt [iweight=tabwt] if (Winner == 1)
		summarize_statorg_capture `JJp1' 4
	}
	else {
		summarize tabwt

		* (1A) Number of Winners
		local col_num = 3
		matrix storedorg_numobs[`JJp1',`col_num'] = r(N)
		matrix storedorg_mean[`JJp1',`col_num'] = 0
		matrix storedorg_variance[`JJp1',`col_num'] = 0
		matrix storedorg_wgtobs[`JJp1',`col_num'] = r(sum)

		* (1B) Avearge Tax Change among Winners
		local col_num = 4
		matrix storedorg_numobs[`JJp1',`col_num'] = r(N)
		matrix storedorg_mean[`JJp1',`col_num'] = 0
		matrix storedorg_variance[`JJp1',`col_num'] = 0
		matrix storedorg_wgtobs[`JJp1',`col_num'] = r(sum)
	}

	count if (Loser == 1)
	if (r(N) > 0) {
		* (1C) Number of Losers
		ratio Loser/ones [iweight=tabwt]
		ratio_statorg_capture `JJp1' 5

		* (1D) Avearge Tax Change among Losers
		*Q* Is this the right statistics for bootstraping since # losers vary across simulations???
		summarize Loser_amt [iweight=tabwt] if (Loser == 1)
		summarize_statorg_capture `JJp1' 6
	}
	else {
		summarize tabwt

		* (1C) Number of Losers
		local col_num = 5
		matrix storedorg_numobs[`JJp1',`col_num'] = r(N)
		matrix storedorg_mean[`JJp1',`col_num'] = 0
		matrix storedorg_variance[`JJp1',`col_num'] = 0
		matrix storedorg_wgtobs[`JJp1',`col_num'] = r(sum)

		* (1D) Avearge Tax Change among Loser
		local col_num = 6
		matrix storedorg_numobs[`JJp1',`col_num'] = r(N)
		matrix storedorg_mean[`JJp1',`col_num'] = 0
		matrix storedorg_variance[`JJp1',`col_num'] = 0
		matrix storedorg_wgtobs[`JJp1',`col_num'] = r(sum)
	}

	* (2) Percent Change in After-Tax Income = Total of diff_AfterTaxIncome / Total of cl_AfterTaxIncome
	ratio diff_AfterTaxIncome/cl_AfterTaxIncome [iweight=tabwt]
	ratio_statorg_capture `JJp1' 7

	* (3) Share of Total Federal Tax Change of group JJ = Total FedTaxChange`JJ' / Total FedTaxChange
	summarize tabwt
	local col_num = 8
	matrix storedorg_numobs[`JJp1',`col_num'] = r(N)
	matrix storedorg_mean[`JJp1',`col_num'] = 1
	matrix storedorg_variance[`JJp1',`col_num'] = 0
	matrix storedorg_wgtobs[`JJp1',`col_num'] = r(sum)

	* (4) Average Federal Tax Change
	summarize FedTaxChange [iweight=tabwt]
	summarize_statorg_capture `JJp1' 9

	* (5) Average Federal Tax Rate Change (% points) = Total FedTaxChange / Total cl_eci_cl
	ratio FedTaxChange/cl_eci_cl [iweight=tabwt]
	ratio_statorg_capture `JJp1' 10

	* (6) Average Federal Tax Rate Under the Proposal = Total alt_burden / Total cl_eci_cl
	ratio alt_burden/cl_eci_cl [iweight=tabwt]
	ratio_statorg_capture `JJp1' 11

matrix list storedorg_numobs
matrix list storedorg_mean
matrix list storedorg_variance
matrix list storedorg_wgtobs

* ----------------------------------
* Collect statistics into one matrix
* ----------------------------------

matrix storedbs_all = storedorg_numobs \ storedorg_mean \ storedorg_variance \ storedorg_wgtobs

* +++++++++++++++++++++++++++++++
* (4) Obtain bootstrap statistics
* +++++++++++++++++++++++++++++++

* Matrix for storing bootstrap estimates
* 12 ECI groups + 1 total = 13 rows, 9 variables

matrix storedbs_numobs   = J(13,1,1), group_num, J(13,9,0) /* number of observations */
matrix storedbs_mean     = J(13,1,2), group_num, J(13,9,0) /* avearge, original */
matrix storedbs_variance = J(13,1,3), group_num, J(13,9,0) /* variance */
matrix storedbs_bsmean   = J(13,1,4), group_num, J(13,9,0) /* avearge, bootstrapped */
matrix storedbs_bias     = J(13,1,5), group_num, J(13,9,0) /* bias */
matrix storedbs_z0       = J(13,1,6), group_num, J(13,9,0) /* z(0) */
matrix storedbs_bsse     = J(13,1,7), group_num, J(13,9,0) /* Bootstrapped standard error */
matrix storedbs_lbcinormal = J(13,1,8), group_num, J(13,9,0) /* Lower bound of the confindence interval -- assuming normal distribution */
matrix storedbs_ubcinormal = J(13,1,9), group_num, J(13,9,0) /* Upper bound of the confindence interval -- assuming normal distribution */
matrix storedbs_lbcip    = J(13,1,10), group_num, J(13,9,0) /* Lower bound of the confindence interval -- raw empirical percentile */
matrix storedbs_ubcip    = J(13,1,11), group_num, J(13,9,0) /* Upper bound of the confindence interval -- raw empirical percentile */
matrix storedbs_lbcibc   = J(13,1,12), group_num, J(13,9,0) /* Lower bound of the confindence interval -- bias-corrected */
matrix storedbs_ubcibc   = J(13,1,13), group_num, J(13,9,0) /* Upper bound of the confindence interval -- bias-corrected */
* Use "ereturn list" to see the list of available statistics after the bootstrap (bs4rw)

foreach JJ of numlist 0/11 {
	local JJp1 = `JJ' + 1

	count if (ECI_lvlgroup == `JJ') & (Winner == 1)
	if (r(N) > 0) {
		* (1A) Number of Winners
		bs4rw, rweights(bsw*): ratio Winner/ones [iweight=tabwt] if (ECI_lvlgroup == `JJ')
		statbs_capture `JJp1' 3

		* (1B) Avearge Tax Change among Winners
		*Q* Is this the right statistics for bootstraping since # winners vary across simulations???
		bs4rw mean=r(mean), rweights(bsw*): summarize Winner_amt [iweight=tabwt] if (ECI_lvlgroup == `JJ') & (Winner == 1)
		statbs_capture `JJp1' 4
	}
	else {
		summarize tabwt if (ECI_lvlgroup == `JJ')

		* (1A) Number of Winners
		local col_num = 3
		matrix storedbs_numobs[`JJp1',`col_num'] = r(N)
		matrix storedbs_mean[`JJp1',`col_num'] = 0
		matrix storedbs_variance[`JJp1',`col_num'] = 0
		matrix storedbs_bsmean[`JJp1',`col_num'] = 0
		matrix storedbs_bias[`JJp1',`col_num'] = 0
		matrix storedbs_z0[`JJp1',`col_num'] = 0
		matrix storedbs_bsse[`JJp1',`col_num'] = 0
		matrix storedbs_lbcinormal[`JJp1',`col_num'] = 0
		matrix storedbs_ubcinormal[`JJp1',`col_num'] = 0
		matrix storedbs_lbcip[`JJp1',`col_num'] = 0
		matrix storedbs_ubcip[`JJp1',`col_num'] = 0
		matrix storedbs_lbcibc[`JJp1',`col_num'] = 0
		matrix storedbs_ubcibc[`JJp1',`col_num'] = 0

		* (1B) Avearge Tax Change among Winners
		local col_num = 4
		matrix storedbs_numobs[`JJp1',`col_num'] = r(N)
		matrix storedbs_mean[`JJp1',`col_num'] = 0
		matrix storedbs_variance[`JJp1',`col_num'] = 0
		matrix storedbs_bsmean[`JJp1',`col_num'] = 0
		matrix storedbs_bias[`JJp1',`col_num'] = 0
		matrix storedbs_z0[`JJp1',`col_num'] = 0
		matrix storedbs_bsse[`JJp1',`col_num'] = 0
		matrix storedbs_lbcinormal[`JJp1',`col_num'] = 0
		matrix storedbs_ubcinormal[`JJp1',`col_num'] = 0
		matrix storedbs_lbcip[`JJp1',`col_num'] = 0
		matrix storedbs_ubcip[`JJp1',`col_num'] = 0
		matrix storedbs_lbcibc[`JJp1',`col_num'] = 0
		matrix storedbs_ubcibc[`JJp1',`col_num'] = 0
	}

	count if (ECI_lvlgroup == `JJ') & (Loser == 1)
	if (r(N) > 0) {
		* (1C) Number of Losers
		bs4rw, rweights(bsw*): ratio Loser/ones [iweight=tabwt] if (ECI_lvlgroup == `JJ')
		statbs_capture `JJp1' 5

		* (1D) Avearge Tax Change among Losers
		*Q* Is this the right statistics for bootstraping since # losers vary across simulations???
		bs4rw mean=r(mean), rweights(bsw*): summarize Loser_amt [iweight=tabwt] if (ECI_lvlgroup == `JJ') & (Loser == 1)
		statbs_capture `JJp1' 6
	}
	else {
		summarize tabwt if (ECI_lvlgroup == `JJ')

		* (1C) Number of Losers
		local col_num = 5
		matrix storedbs_numobs[`JJp1',`col_num'] = r(N)
		matrix storedbs_mean[`JJp1',`col_num'] = 0
		matrix storedbs_variance[`JJp1',`col_num'] = 0
		matrix storedbs_bsmean[`JJp1',`col_num'] = 0
		matrix storedbs_bias[`JJp1',`col_num'] = 0
		matrix storedbs_z0[`JJp1',`col_num'] = 0
		matrix storedbs_bsse[`JJp1',`col_num'] = 0
		matrix storedbs_lbcinormal[`JJp1',`col_num'] = 0
		matrix storedbs_ubcinormal[`JJp1',`col_num'] = 0
		matrix storedbs_lbcip[`JJp1',`col_num'] = 0
		matrix storedbs_ubcip[`JJp1',`col_num'] = 0
		matrix storedbs_lbcibc[`JJp1',`col_num'] = 0
		matrix storedbs_ubcibc[`JJp1',`col_num'] = 0

		* (1D) Avearge Tax Change among Loser
		local col_num = 6
		matrix storedbs_numobs[`JJp1',`col_num'] = r(N)
		matrix storedbs_mean[`JJp1',`col_num'] = 0
		matrix storedbs_variance[`JJp1',`col_num'] = 0
		matrix storedbs_bsmean[`JJp1',`col_num'] = 0
		matrix storedbs_bias[`JJp1',`col_num'] = 0
		matrix storedbs_z0[`JJp1',`col_num'] = 0
		matrix storedbs_bsse[`JJp1',`col_num'] = 0
		matrix storedbs_lbcinormal[`JJp1',`col_num'] = 0
		matrix storedbs_ubcinormal[`JJp1',`col_num'] = 0
		matrix storedbs_lbcip[`JJp1',`col_num'] = 0
		matrix storedbs_ubcip[`JJp1',`col_num'] = 0
		matrix storedbs_lbcibc[`JJp1',`col_num'] = 0
		matrix storedbs_ubcibc[`JJp1',`col_num'] = 0
	}

	* (2) Percent Change in After-Tax Income = Total of diff_AfterTaxIncome / Total of cl_AfterTaxIncome
	bs4rw, rweights(bsw*): ratio diff_AfterTaxIncome/cl_AfterTaxIncome [iweight=tabwt] if (ECI_lvlgroup == `JJ')
	statbs_capture `JJp1' 7

	* (3) Share of Total Federal Tax Change of group JJ = Total FedTaxChange`JJ' / Total FedTaxChange
	*T* For a revenue neutral proposal, FedTaxChange is approximately $0 and thus this ratio is not well-defined!
	bs4rw, rweights(bsw*): ratio FedTaxChange`JJ'/FedTaxChange [iweight=tabwt]
	statbs_capture `JJp1' 8

	* (4) Average Federal Tax Change
	bs4rw mean=r(mean), rweights(bsw*): summarize FedTaxChange [iweight=tabwt] if (ECI_lvlgroup == `JJ')
	statbs_capture `JJp1' 9

	* (5) Average Federal Tax Rate Change (% points) = Total FedTaxChange / Total cl_eci_cl
	bs4rw, rweights(bsw*): ratio FedTaxChange/cl_eci_cl [iweight=tabwt] if (ECI_lvlgroup == `JJ')
	statbs_capture `JJp1' 10

	* (6) Average Federal Tax Rate Under the Proposal = Total alt_burden / Total cl_eci_cl
	bs4rw, rweights(bsw*): ratio alt_burden/cl_eci_cl [iweight=tabwt] if (ECI_lvlgroup == `JJ')
	statbs_capture `JJp1' 11
}

* -------------
* All tax units
* -------------

	local JJp1 = `JJp1' + 1

	count if (Winner == 1)
	if (r(N) > 0) {
		* (1A) Number of Winners
		bs4rw, rweights(bsw*): ratio Winner/ones [iweight=tabwt]
		statbs_capture `JJp1' 3

		* (1B) Avearge Tax Change among Winners
		*Q* Is this the right statistics for bootstraping since # winners vary across simulations???
		bs4rw mean=r(mean), rweights(bsw*): summarize Winner_amt [iweight=tabwt] if (Winner == 1)
		statbs_capture `JJp1' 4
	}
	else {
		summarize tabwt

		* (1A) Number of Winners
		local col_num = 3
		matrix storedbs_numobs[`JJp1',`col_num'] = r(N)
		matrix storedbs_mean[`JJp1',`col_num'] = 0
		matrix storedbs_variance[`JJp1',`col_num'] = 0
		matrix storedbs_bsmean[`JJp1',`col_num'] = 0
		matrix storedbs_bias[`JJp1',`col_num'] = 0
		matrix storedbs_z0[`JJp1',`col_num'] = 0
		matrix storedbs_bsse[`JJp1',`col_num'] = 0
		matrix storedbs_lbcinormal[`JJp1',`col_num'] = 0
		matrix storedbs_ubcinormal[`JJp1',`col_num'] = 0
		matrix storedbs_lbcip[`JJp1',`col_num'] = 0
		matrix storedbs_ubcip[`JJp1',`col_num'] = 0
		matrix storedbs_lbcibc[`JJp1',`col_num'] = 0
		matrix storedbs_ubcibc[`JJp1',`col_num'] = 0

		* (1B) Avearge Tax Change among Winners
		local col_num = 4
		matrix storedbs_numobs[`JJp1',`col_num'] = r(N)
		matrix storedbs_mean[`JJp1',`col_num'] = 0
		matrix storedbs_variance[`JJp1',`col_num'] = 0
		matrix storedbs_bsmean[`JJp1',`col_num'] = 0
		matrix storedbs_bias[`JJp1',`col_num'] = 0
		matrix storedbs_z0[`JJp1',`col_num'] = 0
		matrix storedbs_bsse[`JJp1',`col_num'] = 0
		matrix storedbs_lbcinormal[`JJp1',`col_num'] = 0
		matrix storedbs_ubcinormal[`JJp1',`col_num'] = 0
		matrix storedbs_lbcip[`JJp1',`col_num'] = 0
		matrix storedbs_ubcip[`JJp1',`col_num'] = 0
		matrix storedbs_lbcibc[`JJp1',`col_num'] = 0
		matrix storedbs_ubcibc[`JJp1',`col_num'] = 0
	}

	count if (Loser == 1)
	if (r(N) > 0) {
		* (1C) Number of Losers
		bs4rw, rweights(bsw*): ratio Loser/ones [iweight=tabwt]
		statbs_capture `JJp1' 5

		* (1D) Avearge Tax Change among Losers
		*Q* Is this the right statistics for bootstraping since # losers vary across simulations???
		bs4rw mean=r(mean), rweights(bsw*): summarize Loser_amt [iweight=tabwt] if (Loser == 1)
		statbs_capture `JJp1' 6
	}
	else {
		summarize tabwt

		* (1C) Number of Losers
		local col_num = 5
		matrix storedbs_numobs[`JJp1',`col_num'] = r(N)
		matrix storedbs_mean[`JJp1',`col_num'] = 0
		matrix storedbs_variance[`JJp1',`col_num'] = 0
		matrix storedbs_bsmean[`JJp1',`col_num'] = 0
		matrix storedbs_bias[`JJp1',`col_num'] = 0
		matrix storedbs_z0[`JJp1',`col_num'] = 0
		matrix storedbs_bsse[`JJp1',`col_num'] = 0
		matrix storedbs_lbcinormal[`JJp1',`col_num'] = 0
		matrix storedbs_ubcinormal[`JJp1',`col_num'] = 0
		matrix storedbs_lbcip[`JJp1',`col_num'] = 0
		matrix storedbs_ubcip[`JJp1',`col_num'] = 0
		matrix storedbs_lbcibc[`JJp1',`col_num'] = 0
		matrix storedbs_ubcibc[`JJp1',`col_num'] = 0

		* (1D) Avearge Tax Change among Loser
		local col_num = 6
		matrix storedbs_numobs[`JJp1',`col_num'] = r(N)
		matrix storedbs_mean[`JJp1',`col_num'] = 0
		matrix storedbs_variance[`JJp1',`col_num'] = 0
		matrix storedbs_bsmean[`JJp1',`col_num'] = 0
		matrix storedbs_bias[`JJp1',`col_num'] = 0
		matrix storedbs_z0[`JJp1',`col_num'] = 0
		matrix storedbs_bsse[`JJp1',`col_num'] = 0
		matrix storedbs_lbcinormal[`JJp1',`col_num'] = 0
		matrix storedbs_ubcinormal[`JJp1',`col_num'] = 0
		matrix storedbs_lbcip[`JJp1',`col_num'] = 0
		matrix storedbs_ubcip[`JJp1',`col_num'] = 0
		matrix storedbs_lbcibc[`JJp1',`col_num'] = 0
		matrix storedbs_ubcibc[`JJp1',`col_num'] = 0
	}

	* (2) Percent Change in After-Tax Income = Total of diff_AfterTaxIncome / Total of cl_AfterTaxIncome
	bs4rw, rweights(bsw*): ratio diff_AfterTaxIncome/cl_AfterTaxIncome [iweight=tabwt]
	statbs_capture `JJp1' 7

	* (3) Share of Total Federal Tax Change of group JJ = Total FedTaxChange`JJ' / Total FedTaxChange
	summarize tabwt
	local col_num = 8
	matrix storedbs_numobs[`JJp1',`col_num'] = r(N)
	matrix storedbs_mean[`JJp1',`col_num'] = 1
	matrix storedbs_variance[`JJp1',`col_num'] = 0
	matrix storedbs_bsmean[`JJp1',`col_num'] = 1
	matrix storedbs_bias[`JJp1',`col_num'] = 0
	matrix storedbs_z0[`JJp1',`col_num'] = 0
	matrix storedbs_bsse[`JJp1',`col_num'] = 0
	matrix storedbs_lbcinormal[`JJp1',`col_num'] = 1
	matrix storedbs_ubcinormal[`JJp1',`col_num'] = 1
	matrix storedbs_lbcip[`JJp1',`col_num'] = 1
	matrix storedbs_ubcip[`JJp1',`col_num'] = 1
	matrix storedbs_lbcibc[`JJp1',`col_num'] = 1
	matrix storedbs_ubcibc[`JJp1',`col_num'] = 1

	* (4) Average Federal Tax Change
	bs4rw mean=r(mean), rweights(bsw*): summarize FedTaxChange [iweight=tabwt]
	statbs_capture `JJp1' 9

	* (5) Average Federal Tax Rate Change (% points) = Total FedTaxChange / Total cl_eci_cl
	bs4rw, rweights(bsw*): ratio FedTaxChange/cl_eci_cl [iweight=tabwt]
	statbs_capture `JJp1' 10

	* (6) Average Federal Tax Rate Under the Proposal = Total alt_burden / Total cl_eci_cl
	bs4rw, rweights(bsw*): ratio alt_burden/cl_eci_cl [iweight=tabwt]
	statbs_capture `JJp1' 11

matrix list storedbs_numobs
matrix list storedbs_mean
matrix list storedbs_variance
matrix list storedbs_bsmean
matrix list storedbs_bias
matrix list storedbs_z0
matrix list storedbs_bsse
matrix list storedbs_lbcinormal
matrix list storedbs_ubcinormal
matrix list storedbs_lbcip
matrix list storedbs_ubcip
matrix list storedbs_lbcibc
matrix list storedbs_ubcibc

* ----------------------------------
* Collect statistics into one matrix
* ----------------------------------

matrix storedbs_all = storedbs_all \ storedbs_numobs \ storedbs_mean \ storedbs_variance \ storedbs_bsmean \ storedbs_bias \ storedbs_z0 \ storedbs_bsse \ storedbs_lbcinormal \ storedbs_ubcinormal \ storedbs_lbcip \ storedbs_ubcip \ storedbs_lbcibc \ storedbs_ubcibc

clear
svmat double storedbs_all, names(var)
rename var1 stat_group
rename var2 ECI_lvlgroup /* 99 indicates all tax units */

* Adjust the scale of the values
foreach JJ in 3 5 7 8 10 11 {
	replace var`JJ' = 100*var`JJ'
}

format %32.10f var*
bysort stat_group ECI_lvlgroup: assert (_N == 1) 
save "$EXTRACT_DIR\bs_dist_${PROP_CUR}${num_reps}.dta", replace

log close
