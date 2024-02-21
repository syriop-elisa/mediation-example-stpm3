/*
Author: Elisavet Syriopoulou
Date: February 21st, 2024


This do file includes an example of Stata code on how to obtain the natural direct (NDE) and indirect (NIE) effects using relative survival.

The quantities of interested are discussed in more detail in the following paper:
	Syriopoulou E, Rutherford MJ, Lambert PC. Understanding disparities in cancer prognosis: An extension of mediation
	analysis to the relative survival framework. Biometrical Journal. 2021; 63: 341–353.  https://doi.org/10.1002/bimj.201900355 

Similar code was publised as part of the paper using stpm2 but here the example used the newest command stpm3.

Data required:
- Simulated data based on the colon cancer data (colonsim.dta. The illustrave example is decribed in the paper. 
- A population lifetable with the expected mortality rates stratified by sex, year, age and deprivation is required for this analysis. The data for producing pomort.dta is available at the Office for National Statistics web-site: 		
	https://www.ons.gov.uk/peoplepopulationandcommunity/healthandsocialcare/conditionsanddiseases/datasets/cancersurvivalsmoothedlifetables 
	
The simulated dataset (colonsim.dta) includes the following variables:
* dep: deprivation status (1 for the least deprived group and 5 for the most deprived group)
* agediag: age at diagnosis
* stage: stage at diagnosis (ordinal variable with 1 for the least advanced stage)
* yeardiag: year of diagnosis
* sex: sex (0 for males and 1 for females)
* t: survival time
* dead: status (0 for alive, 1 for dead)
* id: patient id
* diagmonth: month of diagnosis
* datediag: date of diagnosis

The following user-written commands are required and can be installed within Stata :
 - stpm3        //Fits flexible parametric survival models (so called Royston-Parmar models)
 - gensplines   //Generate various types of spline basis functions
 - standsurv    //Postestimation command that calculates various standardized (marginal) measures after fitting a parametric survival model. 
 - erepost      //Module to repost the estimation results, neded for the parametric bootstrap

Run on Stata 18

*/

//Load simulated colon cancer data
use colonsim.dta, clear

//Recode exposure of interest (deprivation) as 1/0
tab dep
gen depind=1 if dep==5 // if most deprived (exposed)
replace depind=0 if dep==1 // if least deprived (unexposed)
drop dep
rename depind dep
tab dep

//Declare survival data
stset t,  failure(dead==1)  id(id) exit(time 3)

//Merge in expected survival using the population lifetable titled popmort.dta	
gen _age = min(int(agediag+_t),99)
gen _year = min(floor(yeardiag + _t),2016)
sort sex  _year  _age dep
merge  m:1 sex  _year  _age dep  using  popmort.dta, keep(match master) keepusing(rate)
sort id


//---------------------------  Net-setting: differences in net probability of death ---------------------------//

******Step 1.  Fit a parametric relative survival model for the time-to event outcome including the exposure (dep), mediator(stage2-stage4), potential confounders (sex, agediag) and appropriate interactions (stage*dep) and time-dependent effects (tvc()).
stpm3 i.dep i.sex @ns(agediag, df(3)) i.stage i.dep##i.stage, scale(lncumhazard) df(4) ///
		tvc(i.sex @ns(agediag, df(3))) dftvc(2) bhazard(rate)
//Store model parameters to use later for the confidence intervals
estimates store surv


//Generate time variable for estimates
range tt 0 3 50

/*
//Standardised survival by stage at diagnosis for the least and most deprived
standsurv , surv timevar(tt) at1(dep 1 stage 1)           ///
            frame(stand_dep_stage, replace)               ///
            at2(dep 1 stage 2)                            ///
            at3(dep 1 stage 3)                            ///
		    at4(dep 1 stage 4)                            ///
            at5(dep 0 stage 1)                            ///
            at6(dep 0 stage 2)                            ///
            at7(dep 0 stage 3)                            ///
			at8(dep 0 stage 4)                            ///
            atvars(surv_1_mostdep surv_2_mostdep surv_3_mostdep surv_4_mostdep ///
					surv_1_leastdep surv_2_leastdep surv_3_leastdep surv_4_leastdep ) 
// Plot the estimates
frame stand_dep_stage {
twoway (line surv_?_leastdep tt, pstyle(p1line..) lpattern(solid dash shortdash)) ///
       (line surv_?_mostdep tt, pstyle(p2line...) lpattern(solid dash shortdash)) , ///
       ytitle("Standardised relative survival") ///
       xtitle("Years since diagnosis") ///
       legend(order(1 "stage I" 2 "stage II" 3 "stage III" 3 "stage IV") rows(4) ring(0) pos(7) ///
       size(small)region(lcolor(white)) symxsize(*0.7) colgap(2.5) rowgap(0.2)) ///
       ylabel(0.0 0.2 0.4 0.6 0.8 1.0,angle(h) format(%9.1f))  ///
       note("Blue lines are high and red lines are low socioecononic position")
}
	
*/

******Step 2. Fit a model for the mediator including the exposure and confounders 
//Generate splines for age (these have also been created from stpm3)
gensplines agediag, df(3) type(ns) gen(ns_age)
local ageknots `r(knots)'

//Fit a multinomial regression model for the low SEP (most deprived/exposed)
quietly: mlogit stage ns_age* i.sex if dep==1
estimates store ph1
//Fit a multinomial regression model for the high SEP (least deprived/unexposed)
quietly: mlogit stage ns_age** i.sex if dep==0
estimates store ph0


******Step 3. For each individual in the study population obtain estimates for the probability of being in a specific level of the mediator at each level of the exposure 
//For low SEP
estimates restore ph1
gen tmpdep = dep
replace dep = 1
predict p11 p12 p13 p14
replace dep = tmpdep
drop tmpdep

//For high SEP
estimates restore ph0
gen tmpdep = dep
replace dep = 0
predict p01 p02 p03 p04
replace dep = tmpdep
drop tmpdep


****** Step 4. Obtain estimates for the standardised relative survival functions at each level of the exposure using the predictions of Step 3 as weights. Contrasts of these estimates can be formed to obtain the NDE and NIE 
	
estimates restore surv

// Total difference in relative survival under everyone being exposed and under everyone being unexposed
standsurv, timevar(tt) failure frame(med, replace)   ///
	         at1(dep 1 stage 1, atindweights(p11))      ///
	         at2(dep 1 stage 2, atindweights(p12))      ///
	         at3(dep 1 stage 3, atindweights(p13))      ///
			 at4(dep 1 stage 4, atindweights(p14))      ///
	         at5(dep 0 stage 1, atindweights(p01))      ///
	         at6(dep 0 stage 2, atindweights(p02))      ///
	         at7(dep 0 stage 3, atindweights(p03))      ///
			 at8(dep 0 stage 4, atindweights(p04))      ///
	         lincom(1 1 1 1 -1 -1 -1 -1) lincomvar(tce_PE) 
		
// NDE		
frame med: capture drop _at*
standsurv, timevar(tt) failure frame(med, merge)   ///
	         at1(dep 1 stage 1, atindweights(p01))      ///
	         at2(dep 1 stage 2, atindweights(p02))      ///
	         at3(dep 1 stage 3, atindweights(p03))      ///
			 at4(dep 1 stage 4, atindweights(p04))      ///
	         at5(dep 0 stage 1, atindweights(p01))      ///
	         at6(dep 0 stage 2, atindweights(p02))      ///
	         at7(dep 0 stage 3, atindweights(p03))      ///
			 at8(dep 0 stage 4, atindweights(p04))      ///
	         lincom(1 1 1 1 -1 -1 -1 -1) lincomvar(tde_PE) 
			 
			 
// NIE
frame med: capture drop _at*
standsurv, timevar(tt) failure  frame(med, merge)   ///
	         at1(dep 1 stage 1, atindweights(p11))      ///
	         at2(dep 1 stage 2, atindweights(p12))      ///
	         at3(dep 1 stage 3, atindweights(p13))      ///
			 at4(dep 1 stage 4, atindweights(p14))      ///
	         at5(dep 1 stage 1, atindweights(p01))      ///
	         at6(dep 1 stage 2, atindweights(p02))      ///
	         at7(dep 1 stage 3, atindweights(p03))      ///
			 at8(dep 1 stage 4, atindweights(p04))      ///
	         lincom(1 1 1 1 -1 -1 -1 -1) lincomvar(tie_PE) 
			 

	
//Make estimates start at 0 for plots
frame med {
  capture drop _at*  
  replace tce_PE=0 if tt==0
  replace tde_PE=0 if tt==0
  replace tie_PE=0 if tt==0
  line tce_PE tde_PE tie_PE tt, sort	
}

//Proportion mediated by stage
frame med: gen prop=tie_PE/tce_PE
frame med: list tt tce_PE tde_PE tie_PE prop if tt==3, noobs

	
frame med {		
	twoway (line tce_PE tde_PE tie_PE tt), ///
		   xtitle("Years since diagnosis", size(medlarge)) ///
		   ytitle("Difference in net probabilities of death", size(medlarge)) ///
		   legend(order(1 "Total" 2 "Natural direct effect" 3 "Natural indirect effect") pos(0) bplacement(neast))  ///
		   name(tce_part_rs, replace)	///
		   ylabel(0.00(0.02)0.10, angle(h) format(%3.2f) )	
}	


****** Obtain confidence intervals
// First write a little program called sampbeta to use later on for drawing the model parameters
// Note: This is general and should work for different types of models.
capture pr drop sampbeta
program sampbeta, 
  syntax name
  tempname b V
  matrix `b' = e(b)
  matrix `V' = e(V)
  local Np = colsof(`b')
  local cnames:  colfullnames `b'
  local rnames:  rowfullnames `b'  
  forvalues i = 1/`Np' {
    tempvar p`i'
    local plist `plist' `p`i''
  }
  tempname m
  frame create `m'
  frame `m' {
    drawnorm `plist', mean(`b') cov(`V') n(1)
    mkmat `plist', matrix(`namelist')
  }
  matrix colnames `namelist' = `cnames'
  matrix rownames `namelist' = `rnames'  
  erepost b = `namelist'  V=`V', noesample 
end

****** Steps 5 and.  For the confidence intervals, repeat from Step 3 for m times while performing parametric bootstrap for the parameter estimates for both models. 

global tcelist
global tielist
global tdelist

// Remember to set a seed to be able to replicate the results
set seed 5693454


// Repreat this m times
// We use 10 iterations here as this is only an example but would need to do far more to get reliable standard errors.
global m 10
forvalues i = 1/$m {

	//Weights for low SEP by drawing the model parameters from a multivariate normal distribution
	preserve
	estimates restore ph1
	sampbeta b
	restore
	capture drop tmpdep
	gen tmpdep = dep
	replace dep = 1
	predict p11_`i'  p12_`i'  p13_`i'  p14_`i' 
	replace dep = tmpdep
	drop tmpdep

	//Weights for high SEP by drawing the model parameters from a multivariate normal distribution
	preserve
	estimates restore ph0
	sampbeta b
	restore
	capture drop tmpdep
	gen tmpdep = dep
	replace dep = 0
	predict p01_`i'  p02_`i'  p03_`i'  p04_`i' 
	replace dep = tmpdep
	drop tmpdep

	//Same as before for the survival model
	preserve
	estimates restore surv
	matrix bsurv = e(b)
	matrix VVVsurv= e(V) 
	sampbeta b
	restore
	
	//Obtain predictions for i bootstrap sample
	
	//Total difference
	standsurv, timevar(tt) failure frame(med, merge)   ///
				 at1(dep 1 stage 1, atindweights(p11_`i'))      ///
				 at2(dep 1 stage 2, atindweights(p12_`i'))      ///
				 at3(dep 1 stage 3, atindweights(p13_`i'))      ///
				 at4(dep 1 stage 4, atindweights(p14_`i'))      ///
				 at5(dep 0 stage 1, atindweights(p01_`i'))      ///
				 at6(dep 0 stage 2, atindweights(p02_`i'))      ///
				 at7(dep 0 stage 3, atindweights(p03_`i'))      ///
				 at8(dep 0 stage 4, atindweights(p04_`i'))      ///
				 lincom(1 1 1 1 -1 -1 -1 -1) lincomvar(tce_`i') 
	 //drop _at predictions
	 frame med: capture drop _at*
    
	//Store estimates
	global tcelist $tcelist tce_`i'

	// NDE
	standsurv, timevar(tt) failure frame(med, merge)   ///
				 at1(dep 1 stage 1, atindweights(p01_`i'))      ///
				 at2(dep 1 stage 2, atindweights(p02_`i'))      ///
				 at3(dep 1 stage 3, atindweights(p03_`i'))      ///
				 at4(dep 1 stage 4, atindweights(p04_`i'))      ///
				 at5(dep 0 stage 1, atindweights(p01_`i'))      ///
				 at6(dep 0 stage 2, atindweights(p02_`i'))      ///
				 at7(dep 0 stage 3, atindweights(p03_`i'))      ///
				 at8(dep 0 stage 4, atindweights(p04_`i'))      ///
				 lincom(1 1 1 1 -1 -1 -1 -1) lincomvar(tde_`i') 
			 		
	//Store estimates
	global tdelist $tdelist tde_`i'		
	frame med: capture drop _at*

	// NIE
	standsurv, timevar(tt) failure  frame(med, merge)   ///
				 at1(dep 1 stage 1, atindweights(p11_`i'))      ///
				 at2(dep 1 stage 2, atindweights(p12_`i'))      ///
				 at3(dep 1 stage 3, atindweights(p13_`i'))      ///
				 at4(dep 1 stage 4, atindweights(p14_`i'))      ///
				 at5(dep 1 stage 1, atindweights(p01_`i'))      ///
				 at6(dep 1 stage 2, atindweights(p02_`i'))      ///
				 at7(dep 1 stage 3, atindweights(p03_`i'))      ///
				 at8(dep 1 stage 4, atindweights(p04_`i'))      ///
				 lincom(1 1 1 1 -1 -1 -1 -1) lincomvar(tie_`i') 
				 
	
	//Store estimates
	global tielist $tielist tie_`i'
	frame med: capture drop _at*


	drop p11_`i' p12_`i' p13_`i' p14_`i' p01_`i' p02_`i' p03_`i' p04_`i'
}
	


	

****** Step 6. Calculate 95% confidence intervalsby using the standard deviation of the estimates obtained from the bootstrap samples.

//For the total difference
frame med {
  egen tce_all_sd=rowsd($tcelist)
  gen tce_all_lci= tce_PE - invnormal(1-(1-c(level)/100)/2) * tce_all_sd
  gen tce_all_uci= tce_PE + invnormal(1-(1-c(level)/100)/2) * tce_all_sd
  //For the NDE
  egen tde_all_sd=rowsd($tdelist)
  gen tde_all_lci= tde_PE - invnormal(1-(1-c(level)/100)/2) * tde_all_sd
  gen tde_all_uci= tde_PE + invnormal(1-(1-c(level)/100)/2) * tde_all_sd
  //For the NIE
  egen tie_all_sd=rowsd($tielist)
  gen tie_all_lci= tie_PE - invnormal(1-(1-c(level)/100)/2) * tie_all_sd
  gen tie_all_uci= tie_PE + invnormal(1-(1-c(level)/100)/2) * tie_all_sd
}


//Make estimates start at 0 for plots
frame med {
  replace tce_PE=0 if tt==0
  replace tce_all_lci=0 if tt==0
  replace tce_all_uci=0 if tt==0
  replace tde_all_lci=0 if tt==0
  replace tde_all_uci=0 if tt==0
  replace tie_all_lci=0 if tt==0
  replace tie_all_uci=0 if tt==0
  replace tde_PE=0 if tt==0
  replace tie_PE=0 if tt==0
}



//Plot with confidence intervals
frame med {

	twoway  (line tce_PE tde_PE tie_PE  tt)                                                   ///
			(rarea tce_all_lci tce_all_uci tt, pstyle(p1line) color(%30))      ///
			(rarea tde_all_lci tde_all_uci tt, pstyle(p2line) color(%30))      ///
			(rarea tie_all_lci tie_all_uci tt, pstyle(p3line) color(%30)),      ///                                                       
		   legend(order(1 "Total" 2 "Natural direct effect" 3 "Natural indirect effect") pos(0) bplacement(neast))  ///
			xtitle("Years since diagnosis", size(medlarge))                    ///
			ytitle("Difference in net probabilities of death", size(medlarge)) ///
			name(all_ci, replace)                                           ///
			ylabel(0.00(0.02)0.10, angle(h) format(%3.2f))	                 
        
}





//---------------------------  All-cause setting: differences in all-cause probability of death using the relative survival framework  ---------------------------//

estimates restore surv

//Need to incorporate the expected mortality rates using option expsurv()

//Total difference in all-cause survival (keeping exposure to the observed level at the expected mortality rates)
frame med: capture drop _at*
standsurv, timevar(tt) failure frame(med, merge)   ///
	         at1(dep 1 stage 1, atindweights(p11))      ///
	         at2(dep 1 stage 2, atindweights(p12))      ///
	         at3(dep 1 stage 3, atindweights(p13))      ///
			 at4(dep 1 stage 4, atindweights(p14))      ///
	         at5(dep 0 stage 1, atindweights(p01))      ///
	         at6(dep 0 stage 2, atindweights(p02))      ///
	         at7(dep 0 stage 3, atindweights(p03))      ///
			 at8(dep 0 stage 4, atindweights(p04))      ///
	         lincom(1 1 1 1 -1 -1 -1 -1) lincomvar(tce_PE_ac) ///
			 expsurv(using(popmort.dta) ///
				datediag(datediag)            ///
				agediag(agediag)        ///
				pmrate(rate)			///
				pmage(_age)				///
				pmyear(_year)            ///
				pmmaxyear(2016)           ///
				pmother(dep sex)	    ///
				at1(.)                              ///  
				at2(.)                              /// 
				at3(.)                              ///  
				at4(.)                              /// 
				at5(.)                              ///  
				at6(.)                              /// 
				at7(.)                              ///  
				at8(.)                              /// 
				)

		
// NDE		
frame med: capture drop _at*
standsurv, timevar(tt) failure frame(med, merge)   ///
	         at1(dep 1 stage 1, atindweights(p01))      ///
	         at2(dep 1 stage 2, atindweights(p02))      ///
	         at3(dep 1 stage 3, atindweights(p03))      ///
			 at4(dep 1 stage 4, atindweights(p04))      ///
	         at5(dep 0 stage 1, atindweights(p01))      ///
	         at6(dep 0 stage 2, atindweights(p02))      ///
	         at7(dep 0 stage 3, atindweights(p03))      ///
			 at8(dep 0 stage 4, atindweights(p04))      ///
	         lincom(1 1 1 1 -1 -1 -1 -1) lincomvar(tde_PE_ac) ///
			 expsurv(using(popmort.dta) ///
				datediag(datediag)            ///
				agediag(agediag)        ///
				pmrate(rate)			///
				pmage(_age)				///
				pmyear(_year)            ///
				pmmaxyear(2016)           ///
				pmother(dep sex)	    ///
				at1(.)                              ///  
				at2(.)                              /// 
				at3(.)                              ///  
				at4(.)                              /// 
				at5(.)                              ///  
				at6(.)                              /// 
				at7(.)                              ///  
				at8(.)                              /// 
				)
			 
			 


// NIE
frame med: capture drop _at*
standsurv, timevar(tt) failure  frame(med, merge)   ///
	         at1(dep 1 stage 1, atindweights(p11))      ///
	         at2(dep 1 stage 2, atindweights(p12))      ///
	         at3(dep 1 stage 3, atindweights(p13))      ///
			 at4(dep 1 stage 4, atindweights(p14))      ///
	         at5(dep 1 stage 1, atindweights(p01))      ///
	         at6(dep 1 stage 2, atindweights(p02))      ///
	         at7(dep 1 stage 3, atindweights(p03))      ///
			 at8(dep 1 stage 4, atindweights(p04))      ///
	         lincom(1 1 1 1 -1 -1 -1 -1) lincomvar(tie_PE_ac) ///
			 expsurv(using(popmort.dta) ///
				datediag(datediag)            ///
				agediag(agediag)        ///
				pmrate(rate)			///
				pmage(_age)				///
				pmyear(_year)            ///
				pmmaxyear(2016)           ///
				pmother(dep sex)	    ///
				at1(.)                              ///  
				at2(.)                              /// 
				at3(.)                              ///  
				at4(.)                              /// 
				at5(.)                              ///  
				at6(.)                              /// 
				at7(.)                              ///  
				at8(.)                              /// 
				)

//Plot
frame med {		
	replace tce_PE_ac=0 if tt==0
	replace tie_PE_ac=0 if tt==0
	replace tde_PE_ac=0 if tt==0
	
	twoway (line tce_PE_ac tde_PE_ac tie_PE_ac tce_PE tde_PE tie_PE tt , lpattern(solid solid solid dash dash dash)), ///
		   xtitle("Years since diagnosis", size(medlarge)) ///
		   ytitle("Difference in all-cause probabilities of death", size(medlarge)) ///
		   legend(order(1 "Total" 2 "Natural direct effect" 3 "Natural indirect effect") pos(0) bplacement(neast))  ///
		   name(med_ac, replace)	///
		   ylabel(0.00(0.02)0.10, angle(h) format(%3.2f) )	
		   
	
}	










//---------------------------  Avoidable deaths  ---------------------------//
// number diagnosed by year
tab dep yeardiag

estimates restore surv


// Total difference
frame med: capture drop _at*
frame change default
standsurv if dep==1, per(2170) timevar(tt) failure frame(med, merge)   ///
	         at1(dep 1 stage 1, atindweights(p11))      ///
	         at2(dep 1 stage 2, atindweights(p12))      ///
	         at3(dep 1 stage 3, atindweights(p13))      ///
			 at4(dep 1 stage 4, atindweights(p14))      ///
	         at5(dep 0 stage 1, atindweights(p01))      ///
	         at6(dep 0 stage 2, atindweights(p02))      ///
	         at7(dep 0 stage 3, atindweights(p03))      ///
			 at8(dep 0 stage 4, atindweights(p04))      ///
	         lincom(1 1 1 1 -1 -1 -1 -1) lincomvar(AD_PE) ///
			 expsurv(using(popmort.dta) ///
				datediag(datediag)            ///
				agediag(agediag)        ///
				pmrate(rate)			///
				pmage(_age)				///
				pmyear(_year)            ///
				pmmaxyear(2016)           ///
				pmother(dep sex)	    ///
				at1(.)                              ///  
				at2(.)                              /// 
				at3(.)                              ///  
				at4(.)                              /// 
				at5(.)                              ///  
				at6(.)                              /// 
				at7(.)                              ///  
				at8(.)                              /// 
				)




// NIE
frame med: capture drop _at*
frame change default
standsurv if dep==1, per(2170) timevar(tt) failure  frame(med, merge)   ///
	         at1(dep 1 stage 1, atindweights(p11))      ///
	         at2(dep 1 stage 2, atindweights(p12))      ///
	         at3(dep 1 stage 3, atindweights(p13))      ///
			 at4(dep 1 stage 4, atindweights(p14))      ///
	         at5(dep 1 stage 1, atindweights(p01))      ///
	         at6(dep 1 stage 2, atindweights(p02))      ///
	         at7(dep 1 stage 3, atindweights(p03))      ///
			 at8(dep 1 stage 4, atindweights(p04))      ///
	         lincom(1 1 1 1 -1 -1 -1 -1) lincomvar(ADs_PE) ///
			 expsurv(using(popmort.dta) ///
				datediag(datediag)            ///
				agediag(agediag)        ///
				pmrate(rate)			///
				pmage(_age)				///
				pmyear(_year)            ///
				pmmaxyear(2016)           ///
				pmother(dep sex)	    ///
				at1(.)                              ///  
				at2(.)                              /// 
				at3(.)                              ///  
				at4(.)                              /// 
				at5(.)                              ///  
				at6(.)                              /// 
				at7(.)                              ///  
				at8(.)                              /// 
				)
			 
frame med {		
	replace AD_PE=0 if tt==0
	replace ADs_PE=0 if tt==0

	
	twoway (line AD_PE ADs_PE tt , lpattern(solid solid solid dash dash dash)), ///
		   xtitle("Years since diagnosis", size(medlarge)) ///
		   ytitle("Number of avoidable deaths") legend(order(1 "Total avoidable deaths" 2 "By stage shifting" ) pos(0) bplacement(neast))  ///
		   name(ad, replace)	///
		   ylabel(0(20)120)	
		   
}	

frame med {	
	li AD_PE ADs_PE tt if tt==3
	
}