global EXTRACT_DIR "K:\TPC\SCRATCH\For Lek\Sloan\data"
global PROGRAM_DIR "K:\TPC\SCRATCH\For Lek\Sloan\Program"

global num_reps = 201

* ========
* Baseline
* ========

do "$PROGRAM_DIR\bs_dist_cl_prog.do"

* ===================
* Tax exempt interest
* ===================

global PROP_CUR "Txexint_2019pr"
do "$PROGRAM_DIR\bs_dist_alt_prog.do"
clear
