*-------------------------------------------------------------------------------
*
* PREPARE CONSUMPTION DATA FROM EVS 2013 DATA
*
* This code uses the 2013 SUF of the EVS (Einkommens- und Verbrauchsstichprobe)
* that is available at the research data center of the German Statistical Office
* https://www.forschungsdatenzentrum.de/de/haushalte/evs
* 
* Last updated: March 7, 2022 by Moritz Kuhn, University of Bonn, 
*               mokuhn@uni-bonn.de
*               Lennard Schlattmann (University of Bonn) provided outstanding 
*               research assistance. 
*-------------------------------------------------------------------------------

************** HOUSEKEEPING ****************************************************

	set more off

	** SET PATHS:
	local wd_mo "C:\Users\mokuhn\sciebo\energy_sources\"
	cd `wd_mo'

	** CONFIGURE LOCAL SETUP:
	local path_in "data/evs2013_aa_gs_hb.dta"
	local path_out "data/consumption13.dta"

	** LOAD DATA:
	use "`path_in'", clear


************* CONSTRUCT GENERAL VARIABLES **************************************

	* WEIGHTS
	gen wgt = EF107
	lab var wgt "Sample weights"
	* QUARTER OF YEAR
	gen quarter = EF6
	lab var quarter "Quarter of year"
	* AGE
	gen age_group = EF40 
	lab var age_group "Age group of main earner"
	gen age_group_large = .
	replace age_group_large = 1 if age_group < 5
	replace age_group_large = 2 if age_group > 4 & age_group < 7
	replace age_group_large = 3 if age_group > 6 & age_group < 10
	replace age_group_large = 4 if age_group > 9 & age_group < 13
	replace age_group_large = 5 if age_group > 12
	
	gen age = 2013 - EF8U3
	cap drop age_group_own
	gen age_group_own = .
	replace age_group_own = 1 if age < 30
	replace age_group_own = 2 if age >= 30 & age < 40
	replace age_group_own = 3 if age >= 40 & age < 45
	replace age_group_own = 4 if age >= 45 & age < 50
	replace age_group_own = 5 if age >= 50 & age < 55
	replace age_group_own = 6 if age >= 55 & age < 60
	replace age_group_own = 7 if age >= 60 & age < 65
	replace age_group_own = 8 if age >= 65 & age < 70
	replace age_group_own = 9 if age >= 70
	
	* HOMEOWNER
	gen renter = 0
	replace renter = 1 if EF20 == 3
	gen home_owner = 0
	replace home_owner = 1 if EF20 == 1 | EF20 == 2
	
	* HOUSEHOLD SIZE:
	gen size = EF7
	lab var size "Household size"
	gen children = .
	replace children = 0 if inlist(EF39,1,2,9,10,21,22)
	replace children = 1 if inlist(EF39,3,11,12,23,24)
	replace children = 2 if inlist(EF39,5,13,14,25,26)
	replace children = 3 if inlist(EF39,7,15,16,27,28)
	replace children = 4 if inlist(EF39,  17,18)
	replace children = 5 if inlist(EF39,  19,20)
	lab var children "Number of children < 18"
	gen size_3 = .
	replace size_3 = min(size,3) 
	
	* CITY-SIZE
	gen city_size = .
	replace city_size = 1 if EF4 == 1 | EF4 == 2 | EF4 == 6
	replace city_size = 2 if EF4 == 3 | EF4 == 8
	replace city_size = 3 if EF4 == 4 | EF4 == 9
	replace city_size = 4 if EF4 == 5
	label define city_size_label 1 "< 20.000" 2 "between 20.000 and 100.000" 3 "between 100.000 and 500.000" 4 "> 500.000"
	label value city_size city_size_label
	* EF == 7, not defined in data set!
	* Problem: Not all states report numbers according to these bins
	*          -> Dresden (525.000) and Leipzig (521.000) are categorized in bin 3. 
	*	       -> Erfurt (203.000), Rostock (203.000), Potsdam (159.000) Jena (107.000) are in bin 2.
	
	* GENDER
	gen main_earner_male = 0
	replace main_earner_male = 1 if EF8U2 == 1
	
	* SOCIAL POSITION
	forvalues i = 1/8 {
		local t = `i' + 7
		gen social_pos_`i' =  EF`t'U8 
		}
	
	* OCCUPATION
	rename EF8U19 occupation
	
	* Generate industry for unemployed
	replace occupation = 23 if EF8U8 == 9 
	* Generate industry for "Rentner"
	replace occupation = 24 if EF8U8 == 10
	* Generate industry for "Pensionär"
	replace occupation = 25 if EF8U8 == 11
	
	* EDUCATION
	gen education = EF8U7
	
	* Highest school degree
	gen education_school = . 
	replace education_school = 1 if EF8U7 < 5 & EF8U7 != .
	replace education_school = 2 if EF8U7 > 4 & EF8U7 != .
	replace education_school = 3 if education > 10 & education != .
	
	* NUMBER OF EMPLOYED PEOPLE
	forvalues i = 1/8 {
		local t = `i' + 7
		gen emp_`i' = 0 if social_pos_`i' != 10 & social_pos_`i' != 11
		replace emp_`i' = 1 if EF`t'U13 > 0 & EF`t'U13 != .
		}
		
	gen number_employed = emp_1 + emp_2 + emp_3 + emp_4 + emp_5 + emp_6 + emp_7 + emp_8
	gen one_earner_HH = 0 if emp_1 != .
	replace one_earner_HH = 1 if number_employed == 1
	
	* Unemployed
	gen unemployed = 0
	replace unemployed = 1 if social_pos_1 == 9
	
	* LIVING SPACE:
	gen living_space = EF21
	lab var living_space "Living Space (in square meters)"
	gen living_space_pc = living_space/size
	lab var living_space_pc "Living Space per HH member (in square meters)"	
		
************ CONSTRUCT INCOME AND EXPENDITURE VARIABLES ************************

	* CONSUMPTION: 
	* Define categories which are not coded as consumption in EVS data (insurance and further consumption)
	* Private Insurances
	gen c_insurance = EF98
	* Further expenditures
	gen c_further = EF530 + EF472 + EF473 + EF474 + EF476 + EF531 // ground rent (Erbpacht)(EF530); membership fees (EF472); money donations (EF473); 
																  // voluntary alimony (EF474); gambling costs (EF476 "Spieleinsätze"); further expenditures (EF531)
	** Total consumption
	gen c_all = EF89 - EF77 + c_insurance + c_further // total private consumption (EF89) - imputed rents (EF77) + additional consumption elements
	lab var c_all "Total consumption minus imputed rents plus insurance and further expenditures (EUR)"
	gen c_all_imp_rents = c_all + EF77
	gen c_all_no_insurance = c_all - c_insurance
	gen c_all_no_ins_no_exp = EF89 - EF77
	
	* Consumption categories
	gen c_food = EF73 + EF74 
	lab var c_food "Food consumption (EUR)"
	gen c_house = EF76 + EF78 + EF79 + EF530 // including rent, maintenance, energy and ground rent (Erbpacht)
	lab var c_house "Housing expenditure (EUR)"
	gen c_h_imp_rents = c_house + EF77
	gen c_h_credit = c_house + EF102
	gen c_cloth = EF75 
	lab var c_cloth "Apparel consumption (EUR)"
	gen c_transp = EF82
	lab var c_transp "Transportation expenditure (EUR)"
	gen c_comm = EF84
	lab var c_comm "Communication expenditure (EUR)"
	gen c_leisure = EF85 + EF87 + EF472 + EF476  // including membership fees and gambling costs (EF476 "Spieleinsätze")
	lab var c_leisure "Leisure consumption (EUR)"
	
	* Residual consumption
	gen c_other              = c_all              - c_food - c_house - c_cloth - c_transp - c_comm - c_leisure - c_insurance
	lab var c_other "Other consumption (health care, household items, education, other goods and services) (EUR)"
	gen c_other_no_insurance = c_all_no_insurance - c_food - c_house - c_cloth - c_transp - c_comm - c_leisure
	gen c_all_credit         = c_all + EF102
	
	* EXPENDITURES:
	gen exp_further_taxes = EF96
	lab var exp_further_taxes "Inheritance and gift taxes, etc."
	gen exp_voluntary_insurance = EF98
	lab var exp_voluntary_insurance "Voluntary contributions for private health insurance" 
	gen exp_further_transfers = EF100
	lab var exp_further_transfers "Membership fees, donations, alimony payments" 
	gen exp_accumulate_wealth = EF101
	lab var exp_accumulate_wealth "Buying/renovating houses, buying shares/stocks/gold etc."
	gen exp_repay_loans = EF102
	lab var exp_repay_loans "Repaying of house/consumption loans and interests"
	gen exp_further = EF103
	lab var exp_further "Further expenditures including exp. for business purposes and ground rent (Erbpacht)"
	gen expenditures_all = c_all + exp_further_taxes + exp_voluntary_insurance + exp_further_transfers
	lab var expenditures_all "Total expenditures (consumption)"
	
	* INCOME:
	* Total income
	gen income_total = EF72
	lab var income_total "Household total income, quarterly (EUR)"
	* Net household income
	gen income_net = EF62
	lab var income_net "Household net income, quarterly (EUR)"
	* Net disposable household income as defined in EVS
	gen income_net_disp_predefined = EF65
	lab var income_net_disp_predefined "Disposable household net income, quarterly (EUR), as predefined in EVS"
	* Own net disposable household income
	gen income_net_disp = EF65 - exp_further_taxes - EF237U1 - EF237U2 - EF237U3 - EF237U4 - EF237U5 - EF237U6  ///
	                     - EF475 - EF477 - EF529 - EF238U1 - EF238U2 - EF238U3 - EF238U4 - EF238U5 - EF238U6
	replace income_net_disp = . if income_net_disp < 1000 // only a few outliers discarded
	lab var income_net_disp "Own disposable household net income including credit repayments, quarterly (EUR)"	
	
*********** COMPUTE QUINTILES OF INCOME ****************************************
	
	* DECILES UNCONDITIONALLY ON NUMBER OF HH MEMBERS
	xtile income_total_no_size = income_total [pw = wgt], nq(5)
	xtile income_net_no_size = income_net [pw = wgt], nq(5)
	xtile income_dec_no_size = income_net_disp [pw = wgt], nq(5)
	
	label define income_quantiles 1 "bottom 20%" 2 "20% - 40%" 3 "40% - 60%" 4 "60% - 80%" 5 "top 20%"	
	label values income_total_no_size income_net_no_size income_dec_no_size income_quantiles
	
	forvalues hhsize = 1(1)3 {
		xtile income_total_`hhsize' = income_total if size_3 == `hhsize' [pw = wgt], nq(5)
		xtile income_net_`hhsize' = income_net if size_3 == `hhsize' [pw = wgt] , nq(5)
		xtile income_dec_`hhsize' = income_net_disp if size_3 == `hhsize' [pw = wgt], nq(5) 
		
		label values income_total_`hhsize' income_net_`hhsize' income_dec_`hhsize' income_quantiles
		}
	

*********** Energy Sources *****************************************************

	drop if EF23 == 1 // drop those heating with electricity
	gen energy_gas        = EF317
	gen energy_oil        = EF320
	gen energy_coal_wood  = EF321
	gen energy_warm_water = EF323	
	gen energy_dist_heat  = EF324
	gen energy_fuel       = EF383
	gen energy_total      = energy_gas + energy_oil + energy_coal_wood + energy_warm_water + energy_dist_heat + energy_fuel

	* Main energy source
	gen energy_main_gas       = 0
	gen energy_main_oil       = 0 
	gen energy_main_coal_wood = 0 
	gen energy_main_other     = 0 
	gen energy_main_no_info   = 0 
	replace energy_main_gas       = 1 if EF23 == 2
	replace energy_main_oil       = 1 if EF23 == 3
	replace energy_main_coal_wood = 1 if EF23 == 4
	replace energy_main_other     = 1 if EF23 == 5
	replace energy_main_no_info   = 1 if EF23 == 0

	* Umlagen
	gen umlagen_dist_heat = EF327
	gen umlagen_gas = EF329
	gen umlagen_oil = EF330
	gen umlagen_other = EF331 + EF332
	gen umlagen_total = umlagen_dist_heat + umlagen_gas + umlagen_oil + umlagen_other
	
	gen expenditure_gas = energy_gas + umlagen_gas
	gen expenditure_oil = energy_oil + umlagen_oil
	gen expenditure_dist_heat = energy_dist_heat + umlagen_dist_heat
	gen expenditure_coal_wood = energy_coal_wood

	/* Add in costs for warm water */
	replace expenditure_gas = expenditure_gas + energy_warm_water if energy_main_gas == 1
	replace expenditure_oil = expenditure_oil + energy_warm_water if energy_main_oil == 1
	replace expenditure_coal_wood = expenditure_coal_wood + energy_warm_water if energy_main_coal_wood == 1
	replace expenditure_dist_heat = expenditure_dist_heat + energy_warm_water if inlist(EF23,0,5)
	
	/* Total energy expenditure */
	gen expenditure_total = expenditure_gas + expenditure_oil + expenditure_coal_wood + expenditure_dist_heat
	gen expenditure_fuel = energy_fuel
	
	/* Share in total expenditure and income */
	foreach C in gas oil coal_wood dist_heat total fuel {
		/* Consumption */
		gen energy_share_c_`C' = expenditure_`C'/c_all * 100
		/* Income */
		gen energy_share_y_`C' = expenditure_`C'/income_net * 100
		quietly : sum energy_share_c_`C' [aw = wgt] 
		local consshare = `r(mean)' 
		quietly : sum energy_share_y_`C' [aw = wgt] 
		local incomeshare = `r(mean)' 
		display %22s "`C' (cons/inc): " %6.1f `consshare' " % " %6.1f `incomeshare' " % "
		}
		

gen heatingtype = EF23
recode heatingtype (2 = 1) (3 = 2) (4 = 3) (5 0 = 4)

label define heatingtype_lbl 1 "Gas" 2 "Oil" 3 "Coal & Wood" 4 "District heating & other" 
label define heatingtype_lbl_short 1 "G" 2 "O" 3 "CW" 4 "D&O" 

tabstat expenditure_* [aw = wgt] , by(income_net_no_size)

local i = 1
foreach shi in "c" "y" { 
	label values heatingtype heatingtype_lbl
	
	if(`i' == 1) {
		local estr = "expenditure"
		}
	else {
		local estr = "net income"
		}
	
	/* All households by type of heating */
	graph bar energy_share_`shi'_total energy_share_`shi'_fuel [aw = wgt] , blabel(total,  format(%4.1f) ) /*
	*/ over(heatingtype) bargap(15) ytitle("`estr' share (in %)") graphregion(color(white)) bgcolor(white) legend(lab(1 "heating") lab(2 "fuel"))  name("F0S`i'", replace)
	graph export "`shi'_household_type_all.png", replace


	label values heatingtype heatingtype_lbl_short

	/* All households by income (fuel and energy) */	 
	graph bar   energy_share_`shi'_gas energy_share_`shi'_oil energy_share_`shi'_coal energy_share_`shi'_dist [aw = wgt], /* blabel(total, format(%4.1f))*/ /*
	*/ over(income_net_no_size) stack legend(lab(1 "Gas") lab(2 "Oil") lab(3 "Coal & Wood") lab(4 "District Heating")) graphregion(color(white)) bgcolor(white) name("F1S`i'", replace)
	graph export "`shi'_household_income_all.png", replace

	/* By heating type and income group */
	graph bar energy_share_`shi'_total if inlist(heatingtype,1,2,4) [aw = wgt] , blabel(total,  format(%4.1f) )/*
	*/ over(heatingtype, gap(80) label(labsize(vsmall)) ) over(income_net_no_size, label(labsize(small))) ytitle("`estr' share (in %)") graphregion(color(white)) bgcolor(white) name("F2S`i'", replace)
	graph export "`shi'_household_income_type_all.png", replace

	forvalues hhsize = 1(1)3 {
		graph bar energy_share_`shi'_total  if size_3 == `hhsize' & inlist(heatingtype,1,2,4) [aw = wgt] , blabel(total,  format(%4.1f) ) /*
	*/ over(heatingtype, gap(80) label(labsize(vsmall)) ) over(income_net_`hhsize', label(labsize(small))) ytitle("`estr' share (in %)") graphregion(color(white)) bgcolor(white) name("F2S`i'H`hhsize'", replace)
		graph export "`shi'_household_income_type_size`hhsize'.png", replace
		}
	local i = `i' + 1
	}