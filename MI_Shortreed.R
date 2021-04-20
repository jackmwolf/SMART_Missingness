##############################################################################
##############################################################################
#####  This code accompanies the following paper as supplementary material:
#####    Shortreed, Laber, Pineau, Stroup, Murphy.
#####	 "A multiple imputation strategy for sequential multiple 
#####	    assignment randomized trials."
#####  It reads in the artificial CATIE data set with missing information
#####   provided as supplementary material and use the imputation 
#####	strategy outlined in the paper to imputed missing information.
#####   This code creates 1 completed artificial CATIE data set with
#####	all non-structurally missing information replaced with draws
#####	from the appropriate posterior predictive distributions
##############################################################################
##### This code was wrtten for R version 3.0.2				 
#####   It requires the mice library, version 2.18,			 
#####	  and the pan library, version 0.9.   				 
##############################################################################
##### The artificial CATIE data was designed using data collected from   
#####   CATIE and contains the following variables in wide format.	 
#####	There are 3 study visits (baseline, follow-up visit 1 [ F-up v1] 
#####   end-of-study/follow-up visit 2 [ES visit])and 2 treatment stages.
#####	While the data were generated based on the CATIE trial the       
#####	artificial data hvar fewer variabels and time points and has     
#####	1500 participants.   	   	     	      	     	 	 
##############################################################################
##############      VARIABLES IN ARTIFICIAL CATIE DATA SET	##############
#####  CATIEID: id variables that indexes each participant		 
#####  EMPLOY: Categorical variable indicating employment status at baseline
#####  	       category values: 'Did Not Work'; 'Did Work'
#####  AGE: Continuous values age at enrollment        
#####  SEX: Gender. 'Female'; 'Male'       
#####  TD: Indicator if participant has diagnosis of Tardive Dyskinsia
#####  	   a movement disorder that is a side effect of many 
#####	   antipsychotic medications. `Yes'; 'No'         
#####  EXACER: Indicator for if the participant had been hospitalized for
#####  	       Schizophrenic symptoms in the year prior to CATIE enrollment  
#####	       'Yes'; 'No'   
#####  SITETYPE: Categorical variable indicating the type of study site
#####  		 individual was enrolled through.
#####		 category values: 'Combo'; 'Private Practice'; 
#####		 'State Mental Health'; 'University Clinic'; 'VA' 
#####  YRS_PRES: Years on prescription medication for schizophrenic symptoms
#####  		 prior to CATIE enrollment   
#####  TREATs1: Initial randomly assigned CATIE treatment
#####  		initial treatment options: 
#####		     `Olanzapine`; `Perphenazine`; `Quetiapine` `Risperidone`  
#####   
#####  BMI.0, BMI.1, BMI.2: Time varying values for body mass index (BMI). 
#####  	      	     	    measured at baseline, F-up v1 and ES visit.      
#####  CS16.0, CS16.1, CS16.2: Time varying categorical variable of illness 
#####  	       	       	      severity, as determing by the study physicial
#####			      category values: 'Not ill'; 'Mildly ill';
#####			      'Moderately ill'; 'Markedly ill'; 'Severely ill' 
#####  EPSMEAN.0, EPSMEAN.1, EPSMEAN.2,: Time-varying continuous value of 
#####  		  movement disorder scale.  Is 0 if particiant does not 
#####		  exhibit any symptomes and the scale value for those who do.  
#####  PANSSTOT.0, PANSSTOT.1, PANSSTOT.2: Time-varying continuous value of
#####  		   PANSS (Positive and Negative Symptom Subscale) scale. This 
#####		   is the primary measurement of schizophrenic symptoms used 
#####		   in clincial practice.   
#####  QOL_TOT.0, END_QOL_TOT.2, QOL_TOT.2: Time varying continuous valued
#####  		  Quality of life scale. This information is collected on
#####		  everyone at baseline and at the end-of-study visits, but 
#####		  is only collected on individuals who chose to transition
#####		  into the second CATIE treatment stage; thus is structurally
#####		  missing on individuals who remain on their initial treatment.
#####  MEDAD03.1, MEDAD03.2: Time-varying CATIE medication adherence 
#####  		  	     Proportion of study pills taken in the past 
#####			     treatment period. Medication  adherence 
#####		  	     at baseline is not defined.  
#####  SWITCH: Binary indicator for if individuals chose to transition in to 
#####  	       second treatment stage at the follow-up visit or remain in 
#####	       initial treatment. 'STAYED';'SWITCHED'     
#####  STAGE.0, STAGE.1: The particiants treatment stage.   
#####  			 Stage.0 is 'Stage 1' for everyone 
#####			 Stage.1 is `Stage 2' for those individuals 
#####			 	 who transitioned into second stage and 
#####			 	 `Stage 1' for those who didn't.   
#####  STAGE2_ARM: For those who switch this is the stage 2 randomization arm.
#####  		   category values: 'Tolerability'; 'Efficacy'
#####  TREATs2: Stage 2 randomized treatment.  This is structurally missing for
#####  		individuals who don't transition into second treatment stage. 
#####  TREAT.1. TREAT.2: CATIE treatments.  TREAT.1 is the initial randomly 
#####  			 assinged stage 1 treatent.  TREAT.2 is the initial  
#####			 CATIE treatment for those who stayed in their first  
#####			 treatment through the end of study and for those  
#####			 who transitioned into stage 2, TREAT.2 is their 
#####  			 Stage 2 randomly assigned treatment. 
##############################################################################
##### The flag Flag_Transition_model_impute determines if transition times
#####   are singly (Flag_Transition_model_impute = FALSE) or multiply
#####   (Flag_Transition_model_impute=TRUE) imputed.  Single
#####   imputations are as described on page x of the paper; multiple
#####   imputations are conducted by imputing a binary `transition indicator'
#####   at each time point.
##############################################################################
#source("CATIE_artificial_impute.r")

# load("Artificial_CATIE_Data.Rdata")

## Required imputation libraries
library(mice)
library(pan)

## This file contains all the functions used in the imputations
source("CATIE_artificial_imputations_utilities.r")


## Flag for multiple (TRUE) or single (FALSE) imputation
Flag_transition_model_impute = TRUE

##############################################################################
####  The remainder of this code implements the imputation strategy for the
####  artificial CATIE data set.This code was adapted from the Shortreed et al.
####  paper.
##############################################################################

### Order variables by amount of missingness within each month;
### put the time independent variables first.
# load("Artificial_CATIE_Data.Rdata") ## loads data.frame named CATIE

#Code is weird if there is no missing baseline info 
#Make one observation in each baseline variable missing so code
#can run

smart.mi = function(dat,mis.bas=FALSE,mis.end=TRUE,mis.switch=TRUE){
  Flag_transition_model_impute = TRUE
  dat$CATIEID = 1:length(dat$PANSSTOT.0)
  CATIE_miss=CATIE=dat
  
if (mis.bas==FALSE){
dat$EPSMEAN.0[1]=NA;dat$BMI.0[2]=NA;dat$CS16.0[3]=NA;dat$QOL_TOT.0[4]=NA
dat$PANSSTOT.0[6]=NA;dat$YRS_PRES[8]=NA}
  
if (mis.end==FALSE){
  dat$EPSMEAN.1[1]=NA;dat$BMI.1[2]=NA;dat$CS16.1[3]=NA;dat$END_QOL_TOT.1[4]=NA
  dat$PANSSTOT.1[6]=NA;dat$MEDAD03.1[7]=NA}
  
if (mis.switch==FALSE){
  dat$SWITCH[c(1,3)]=NA;dat$STAGE2_ARM[c(1,3)]=NA;dat$TREATs2[c(1,3)]=NA;dat$STAGE.1[c(1,3)]=NA}

id_var = "CATIEID"
visit_ids = 0:2
variablesNoMiss = unique(c("AGE","SEX","TD","EXACER", "SITETYPE","TREATs1"))
Base_var_miss = c("YRS_PRES", "BMI","CS16","QOL_TOT","EPSMEAN","PANSSTOT")
FVis_1_lvars = c("BMI","CS16","EPSMEAN","MEDAD03","PANSSTOT")
switch_vars = c("SWITCH","STAGE2_ARM","TREATs2","TREAT.1","STAGE.1")
End_of_stage_vars =  c("BMI","CS16","QOL_TOT","PANSSTOT","EPSMEAN","MEDAD03")
FVis_2_lvars = c("BMI","CS16","QOL_TOT","EPSMEAN","MEDAD03","PANSSTOT")

##############################################################################
####  Reorder variables so data is in time-ordered data structure         
##############################################################################
for(ii in visit_ids){
  temp_now = regexpr(paste(".",ii,sep=""), fixed=T,names(CATIE))
  vis_vars = names(CATIE)[which(temp_now > 0)]
  other_vars = names(CATIE)[-which(temp_now > 0)]
  sum_na = colSums(is.na(CATIE[,vis_vars]))
  iorder = sort(sum_na,index.return=T)
  vis_vars = vis_vars[iorder$ix]
  if(ii > 0){
    temp = vis_vars[-which(vis_vars==paste("PANSSTOT",ii,sep="."))]
    vis_vars = c(temp,paste("PANSSTOT",ii,sep="."))
  }
  CATIE = CATIE[,c(other_vars, vis_vars)]
}

##############################################################################
#### Impute variables at baseline
#### Impute PANSS with a MICE at the baseline visit, as no longitudinal data
##############################################################################
#if (mis.bas==TRUE){
  temp_now = regexpr(paste(".",0,sep=""), fixed=T,names(CATIE))>0
  base_vars_miss = names(CATIE)[which(temp_now > 0)]
  base_vars_miss = base_vars_miss[colSums(is.na(CATIE[,base_vars_miss]))!=0]
  variables2impute = c("YRS_PRES",base_vars_miss)
  
  CATIE = impute_CATIE_variables(CATIE=CATIE, seed=round(runif(1)*32000),
                                 cur_visit_id=0, variables2impute=variables2impute,
                                 id_var=id_var, variablesNoMiss=variablesNoMiss, END=FALSE )
  
  variablesNoMiss = unique(c(variablesNoMiss,variables2impute)) 
# }else{
#     variablesNoMiss =c("AGE","SEX", "TD","EXACER","SITETYPE","TREATs1","YRS_PRES","CS16.0","PANSSTOT.0",
#                        "QOL_TOT.0","EPSMEAN.0","BMI.0")}
##############################################################################
#### Impute Stage transition variables for CATIE participants 
#### If Flag_transition_model_impute = TRUE, then use an imputation model
####   to impute if an individual transitioned
#### If Flag_transition_model_impute = FALSE, impute everyone who
####   dropped out of CATIE to transition immediately into the
####   next treatment stage. See the paper for details.
##############################################################################
CATIE$STAGE.1 = "Stage 1"
if( Flag_transition_model_impute ){
  CATIE = impute_CATIE_stage_transition(CATIE,seed=round(runif(1)*32000),
                                        id_var=id_var, variablesNoMiss=variablesNoMiss )
}else{
  CATIE$SWITCH[is.na(CATIE$SWITCH)] = "SWITCHED"
}
CATIE$STAGE.1[CATIE$SWITCH=="SWITCHED"] = "Stage 2"

##############################################################################
#### Create CATIE end-of-stage variables, for those who transition
#### into stage 2 or who were imputed to transition into stage 2
##############################################################################
for( ii in End_of_stage_vars ){
  if( length(which(names(CATIE)==paste("END_",ii,".",1,sep="")))==0 ){
    # Create end of stage variable if need be
    CATIE[,paste("END_",ii,".",1,sep="")] = CATIE[,paste(ii,1,sep=".")]
    # Make sure that those who stayed in stage 1 have missing
    # end of stage variables
    CATIE[CATIE$SWITCH == "STAYED",paste("END_",ii,".",1,sep="")] = NA
    # Make those who transitioned into stage 2 missing for regularly 
    # scheduled variables, so that regularly scheduled and end-of-stage
    # variables are imputed separately
    CATIE[CATIE$SWITCH == "SWITCHED",paste(ii,1,sep=".")] = NA
  }
}
##############################################################################
#### Impute missing end-of-stage variables at follow-up visit 1
#### - Except PANSS
##############################################################################
## Impute stage 2 randomization choice
variables2impute = "STAGE2_ARM"
CATIE = impute_CATIE_variables(CATIE=CATIE, cur_visit_id=1, id_var=id_var,
                               seed=round(runif(1)*32000), END=FALSE, 
                               variables2impute=variables2impute, 
                               variablesNoMiss=variablesNoMiss )
CATIE$STAGE2_ARM[CATIE$SWITCH=="STAYED"] = NA

## Separate out data set for end of stage variables
CATIE_tol = CATIE[!is.na(CATIE$STAGE2_ARM) & 
                    CATIE$STAGE2_ARM == "Tolerability",]
CATIE_eff = CATIE[!is.na(CATIE$STAGE2_ARM) & CATIE$STAGE2_ARM == "Efficacy",]
CATIE_stay = CATIE[is.na(CATIE$STAGE2_ARM),]
##############################################################################
####  This is a hack to nest predictors for the end-of-stage imputation models
####   within stage 2 randomization arm.  mice passive argument does not work
####   For the actual CATIE imputations as described in the manuscript
####   this code accompanies (see top of code for this information), we 
####   wrote over specific Mice functions to allow interactions in the
####   imputation models.
##############################################################################
##  Impute missing end of stage information for those in the tolerability arm
##   and the efficacy arm separately
variables2impute=c( paste( End_of_stage_vars[
  -which(End_of_stage_vars=="PANSSTOT")],1,sep="."))
variables2impute = c(paste("END_",variables2impute,sep=""))
if( sum(is.na(CATIE_tol[,variables2impute])) > 0){
  CATIE_tol = impute_CATIE_variables(CATIE=CATIE_tol, cur_visit_id=1, END=TRUE,
                                     id_var=id_var, seed=round(runif(1)*32000), 
                                     variables2impute=variables2impute, 
                                     variablesNoMiss=variablesNoMiss )
}
if( sum(is.na(CATIE_eff[,variables2impute])) > 0){
  CATIE_eff = impute_CATIE_variables(CATIE=CATIE_eff, cur_visit_id=1, END=TRUE, 
                                     id_var=id_var, seed=round(runif(1)*32000), 
                                     variables2impute=variables2impute, 
                                     variablesNoMiss=variablesNoMiss )
}
## Create binary symptom variable in data set for those who stay
CATIE_stay = make_binary_side_effect_CATIE( CATIE=CATIE_stay, cur_visit_id=1,
                                            sym_var="END_EPSMEAN",bsym_var="END_bEPSMEAN" )
CATIE_stay = CATIE_stay[,names(CATIE_tol)]
CATIE_eff = CATIE_eff[,names(CATIE_tol)]
CATIE = rbind(CATIE_stay,CATIE_eff,CATIE_tol)
CATIE = CATIE[sort(CATIE$CATIEID,index.return=TRUE)$ix,]

##############################################################################
#### Impute missing scheduled variables at follow-up visit 1
####  - Except PANSS
##############################################################################
variables2impute = c(paste(FVis_1_lvars[-which(FVis_1_lvars=="PANSSTOT")],
                           1,sep="."))
CATIE = impute_CATIE_variables(CATIE=CATIE, seed=round(runif(1)*32000),
                               cur_visit_id=1, id_var=id_var, END=FALSE , 
                               variables2impute=variables2impute, 
                               variablesNoMiss=variablesNoMiss)
variablesNoMiss = unique(c(variablesNoMiss,variables2impute))

##############################################################################
#### Combine end of stage & regularly schedule variable for PANSS imputation
##############################################################################
for(ii in c("BMI","CS16","bEPSMEAN","EPSMEAN","MEDAD03","QOL_TOT","PANSSTOT")){
  CATIE[CATIE$SWITCH == "SWITCHED",paste(ii,1,sep=".")] = 
    CATIE[CATIE$SWITCH == "SWITCHED",paste("END_",ii,".",1,sep="")]
}
# drop end of stage variables
CATIE = CATIE[, regexpr("END_",names(CATIE)) < 0 ]
CATIE$QOL_TOT.1[CATIE$SWITCH == "STAYED"] = 0
variablesNoMiss = c(variablesNoMiss, "QOL_TOT.1")
variablesNoMiss = variablesNoMiss[regexpr("END_",variablesNoMiss) < 0]

##############################################################################
#### Impute PANSS at visit 1, both regularly scheduled & end of phase variables
##############################################################################
time_varying_PANSS_predictors = c("BMI","CS16","QOL_TOT","EPSMEAN","bEPSMEAN",
                                  "MEDAD03","TREAT")
CATIE = CATIE_impute_PANSS(CATIE, cur_visit_id=1, id_var=id_var, END=TRUE, 
                           SCH=TRUE,variablesNoMiss=variablesNoMiss,
                           time_varying_PANSS_predictors=
                             time_varying_PANSS_predictors)
variablesNoMiss = c(variablesNoMiss,"PANSSTOT.1")

##############################################################################
####  Impute treatment for those individuals who transitioned into stage 2
##############################################################################
CATIE = CATIE_impute_stage2_treatment(CATIE)

##############################################################################
#### Impute missing scheduled variables at follow-up visit 2
####   - Except PANSS
##############################################################################
variables2impute = c(paste(FVis_2_lvars[-which(FVis_2_lvars=="PANSSTOT")],
                           2,sep="."))
CATIE = impute_CATIE_variables(CATIE=CATIE, seed=round(runif(1)*32000),
                               cur_visit_id=2, id_var=id_var, END=FALSE, 
                               variables2impute=variables2impute, 
                               variablesNoMiss=variablesNoMiss )
variablesNoMiss = unique(c(variablesNoMiss,variables2impute))

##############################################################################
#### Impute PANSS at visit 2, only regularly scheduled in this visit 
##############################################################################
time_varying_PANSS_predictors = c("BMI","CS16","EPSMEAN","bEPSMEAN",
                                  "MEDAD03","TREAT")
CATIE = CATIE_impute_PANSS(CATIE, cur_visit_id=2, id_var=id_var, END=FALSE, 
                           SCH=TRUE,time_varying_PANSS_predictors=
                             time_varying_PANSS_predictors,
                           variablesNoMiss=variablesNoMiss )

## Remove "bEPSMEAN" variables,  they were just used for imputations.
temp_remove = which(regexpr("bEPSMEAN",names(CATIE)) > 0)
CATIE = CATIE[,-temp_remove]

## Replace QOL_TOT variable with NA for those who did not switch
CATIE$END_QOL_TOT.1 = CATIE$QOL_TOT.1
CATIE$END_QOL_TOT.1[CATIE$SWITCH=="STAYED"] = NA
CATIE = CATIE[,-which(names(CATIE)=="QOL_TOT.1")]

## Recode imputed treatments and stages as factors
CATIE$TREAT.1 = factor(CATIE$TREAT.1)
CATIE$TREAT.2 = factor(CATIE$TREAT.2)
CATIE$STAGE.0 = factor(CATIE$STAGE.0)
CATIE$STAGE.1 = factor(CATIE$STAGE.1)

## Reorder CATIE by time.
CATIE_order = c("CATIEID", "AGE", "EMPLOY", "EXACER", "SEX","SITETYPE", 
                "TD", "YRS_PRES", "TREATs1", "BMI.0", "CS16.0", "EPSMEAN.0",
                "PANSSTOT.0", "QOL_TOT.0", "STAGE.0", "TREAT.1",
                "SWITCH", "STAGE.1", "TREATs2" , "STAGE2_ARM", "TREAT.2", 
                "BMI.1", "CS16.1", "EPSMEAN.1", "MEDAD03.1", 
                "PANSSTOT.1", "END_QOL_TOT.1", "BMI.2", "CS16.2", "EPSMEAN.2",  
                "MEDAD03.2", "PANSSTOT.2", "QOL_TOT.2")
setdiff(names(CATIE),CATIE_order)
setdiff(CATIE_order, names(CATIE))
CATIE = CATIE[,CATIE_order]
return(CATIE)
}

#  ## Assess the imputations by plotting the observed versus the imputed values
#  pdf("CATIE_imputations_check.pdf")
#  qqplot(CATIE[is.na(CATIE_miss$YRS_PRES),"YRS_PRES"],
# 	CATIE[!is.na(CATIE_miss$YRS_PRES),"YRS_PRES"],
# 	main="YRS_PRES",xlab="Imputed values",
# 	ylab="Observed values")
#  for(iv in End_of_stage_vars){
#   par(mfrow = c(1,3))
#   for(ii in 0:2){
#    varn = paste(iv,ii,sep=".")
#    if( sum(names(CATIE)==varn) | (varn == "QOL_TOT.1") ){
#     if( varn == "QOL_TOT.1" ){
#      varn = paste("END",varn,sep="_")
#     }  
#     qqplot(CATIE[is.na(CATIE_miss[,varn]),varn],
# 	CATIE[!is.na(CATIE_miss[,varn]),varn],
# 	main=varn,xlab="Imputed values",
# 	ylab="Observed values")
#    }
#   }
#  }
#  dev.off()
