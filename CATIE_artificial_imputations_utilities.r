############################################################################
### List of all functions contained in this file:
###  Details of what the functions do are with the function definition
############################################################################
## impute_CATIE_variables=function( CATIE, seed, cur_visit_id, print_flag=F, 
##				 id_var=id_var, variables2impute, 
##				 variablesNoMiss, END=FALSE )
## define_impute_formulas = function(CATIE_imp, vars_NA, p, variablesNoMiss)
## make_binary_side_effect_CATIE = function(CATIE, cur_visit_id,
##			      sym_var="EPSMEAN",bsym_var="bEPSMEAN")
## clean_up_side_effects_CATIE= function(CATIE, cur_visit_id, END){
## CATIE_get_bounds = function( CATIE, vars_list )
## CATIE_impute_stage2_treatment=function( CATIE )
## CATIE_impute_PANSS=function( CATIE, cur_visit_id, id_var, END=TRUE,
##			    time_varying_PANSS_predictors, SCH=TRUE,  
##  	  	            variablesNoMiss, seed )
## create_predictors4pan = function(CATIE_local,id_var,predictorsNoNest,
##		      	              END=END,SCH=SCH,predictorsNest)
##
## CATIE_impute_stage_transition = function(CATIE,variablesNoMiss,seed
###########################################################################

###########################################################################
## This function imputes CATIE variables: variables2impute
## 	using variablesNoMiss as predcitors.
###########################################################################
impute_CATIE_variables=function( CATIE, seed, cur_visit_id, print_flag=F, 
				 id_var=id_var, variables2impute, 
				 variablesNoMiss, END=FALSE ){

 ## Creates a binary variable for side effect that indicates if the symptom is
 ## present or not.  Side-effect variable is semi-continuous variable, 
 ## so break it up into binary and contiuous portion for iputations
 if( sum(variables2impute == paste("END_EPSMEAN",cur_visit_id,sep=".") ) != 0){
  CATIE = make_binary_side_effect_CATIE( CATIE=CATIE, cur_visit_id,
		      sym_var="END_EPSMEAN",bsym_var="END_bEPSMEAN" )
  		      variables2impute = c(variables2impute, 
		      paste("END_bEPSMEAN",cur_visit_id,sep="."))
 } 
 if( sum(variables2impute == paste("EPSMEAN",cur_visit_id,sep=".") ) != 0){
  CATIE = make_binary_side_effect_CATIE( CATIE=CATIE, cur_visit_id,
		      sym_var="EPSMEAN",bsym_var="bEPSMEAN" )
  variables2impute = c(variables2impute, 
  		   paste("bEPSMEAN",cur_visit_id,sep="."))
 } 

 ## Separate CATIE into data set involved in imputations, variables that either
 ##   need to be imputed or variables  used as predictors for other imputaions
 ## And a data set of variables not involved in the imputations at this stage
 ##  all variables that were collected after cur_visit_id
 ### if need to imputed EPSMEAN, make sure that bEPSMEAN is added to the list
 ### to impute as well.
 CATIE_imp = CATIE[,unique(c(variablesNoMiss,variables2impute))]
 CATIE_other = CATIE[,unique( c(id_var,
  	         setdiff(names(CATIE),c(variablesNoMiss,variables2impute)) ))]

 ## Make predictor matrix so only variables collected prior current visit are 
 ## used as predictors
 p = matrix(0,nrow=length(names(CATIE_imp)),ncol=length(names(CATIE_imp)))
 colnames(p) = rownames(p) = names(CATIE_imp)
 p[variables2impute,variablesNoMiss] = 1
 diag(p) = 0
 
 ## Create formulas for mice to impute missing values
 temp_list = define_impute_formulas( CATIE_imp, vars_NA=variables2impute, p, 
  		    	      variablesNoMiss=c(id_var,variablesNoMiss) ) 
 impute_formulas = temp_list$impute_formulas 
 if(length(which(names(impute_formulas)==id_var)) > 0){
  impute_formulas = impute_formulas[-which(names(impute_formulas)==id_var)] 
 }
 CATIE_imp = temp_list$CATIE_imp
 visit_order = 1:ncol(CATIE_imp)
 visit_order = visit_order[colSums(is.na(CATIE_imp))>0]
 for(ii in names(CATIE_imp)){
  if(is.factor(CATIE_imp[,ii])){
   CATIE_imp[,ii] = factor(CATIE_imp[,ii])
  }
  if( sum(!is.na(CATIE_imp[,ii])) == 1){
   CATIE_imp[is.na(CATIE_imp[,ii]),ii] = CATIE_imp[!is.na(CATIE_imp[,ii]),ii]
  }
 }

 ### Use MICE to impute missing scheduled variabels at this visit
 tmi = mice(data=CATIE_imp, m = 1, method=impute_formulas, 
       predictorMatrix=p, visitSequence=visit_order, printFlag=F, 
       diagnostics=F, seed=seed, 
       defaultMethod=c("norm","logreg","polyreg"))
 CATIE_imp = complete(tmi)

 ### add the id variable back in so CATIE can be merged back in
 CATIE_imp = cbind(CATIE[,id_var],CATIE_imp)
 names(CATIE_imp)[1] = "CATIEID"
 CATIE = merge(CATIE_other, CATIE_imp)
 if(sum(regexpr("EPSMEAN",variables2impute)>0)>0){
  CATIE = clean_up_side_effects_CATIE( CATIE=CATIE, cur_visit_id, END=END )
 }
 ##################################################################
 ### Make sure that all the variables are within the possible bounds
 ##################################################################
 temp = CATIE_get_bounds(vars_list=c(variablesNoMiss,variables2impute))
 upper_bounds = temp[["upper_bounds"]]  
 lower_bounds = temp[["lower_bounds"]]  
 for( ii in names(upper_bounds[!is.na(upper_bounds)]) ){
  if( sum(CATIE[,ii] > upper_bounds[ii]) > 0 ){
   CATIE[CATIE[,ii] > upper_bounds[ii],ii] = upper_bounds[ii]
  }
 }
 for(ii in names(lower_bounds[!is.na(lower_bounds)]) ){
  if( sum(CATIE[,ii] < lower_bounds[ii]) > 0 ){
   CATIE[CATIE[,ii] < lower_bounds[ii],ii] = lower_bounds[ii]
  }
 }
 return(CATIE)
}

##############################################################################
##### The function returns a vector containing the imputation method     #####
#####   for all variables in the CATIE data set.  For variables 	 #####
#####   involving interactions with variables, the vector entry 	 #####
#####	contains a string with the formula for the prediction		 #####
#####   For variables not involving interactions the element contains    #####
#####	"norm" for continuous variables, "logreg" for binary variables   #####
#####	and "polyreg" for categorical variables with more than 2 	 #####
#####	categories    	  	      		     	       		 #####
##############################################################################
#####  CATIE: a CATIE data set as input to CATIE_impute 		 #####
#####  vars_NA: a vector of length ncol(CATIE) indicating the number 	 #####
#####  		of missing values present for each variable   		 #####
#####  p: the prediction matrix to be input into mice			 #####
##############################################################################
define_impute_formulas = function(CATIE_imp, vars_NA, p, variablesNoMiss){

 temp_NA = names(CATIE_imp)[colSums(is.na(CATIE_imp))>0]
 vars_NA = unique(c(vars_NA,temp_NA))
 impute_formulas = rep("",dim(CATIE_imp)[2])
 names(impute_formulas) = names(CATIE_imp)

 # for each of the variables that need to be imputed loop throw and make the 
 # forumla to input input MICE
 for( i in vars_NA ){
  this_visit_id = ""
  ## get the visit number of this variable
  if( substr( i, start=nchar(i)-3, stop=nchar(i)-3) == "."){
    this_visit_id = as.numeric(substr(i,start=nchar(i)-2,stop=nchar(i)))
  }
  if( substr( i, start=nchar(i)-4, stop=nchar(i)-4) == "."  ){
    this_visit_id = as.numeric(substr(i,start=nchar(i)-3,stop=nchar(i)))
  }

  if( (sum(is.na(CATIE_imp[,i]))==0) & (sum(vars_NA==i)>0) ){
   ### if no variables with missing data print warning and leave function
   print(paste("Variable:",i,"something wrong here, nothing to impute..."))
  }else if( (sum(is.na(CATIE_imp[,i]))==0) ){
   ### if there is nothing to impute, formula is be ""
   impute_formulas[i] = ""
  }else if( length(unique(CATIE_imp[,i]))==2 ){
   ### This means there is only one uniue non missing value
   ### fill in all the missing values with this one value and move on...
   t = unique(CATIE_imp[,i])[!is.na(unique(CATIE_imp[,i]))]
   CATIE_imp[is.na(CATIE_imp[,i]),i] = t
   impute_formulas[i] = ""
  } else {
   # Use norm for continuous variables
   if(!is.factor(CATIE_imp[,i])){  
      impute_formulas[i] = "norm" 
   }else{
    # Use logreg or polyreg as appropriate
    CATIE_imp[,i] = factor(CATIE_imp[,i])
    if( nlevels(CATIE_imp[,i]) == 2 ){  
     impute_formulas[i] = "logreg"
    }else{  
     impute_formulas[i] = "polyreg"
    }
   }
  }# close missing values present loop
  #  print(impute_formulas[i])
 }# close for loop over variables
 return(list(impute_formulas=impute_formulas,CATIE_imp=CATIE_imp))
} # close make impute formulas function loop

###########################################################################
## The movement disorder scores in CATIE have a semi-continuous distribution
##  Most people have no movement disorder, i.e. a value of 0 on the score
##  Make a binary variable that indicated if an individual has a movement
##  disorder or not in addition to the continuous valued score for those
##  who have a movement disorder
###########################################################################
make_binary_side_effect_CATIE = function(CATIE, cur_visit_id,
			      sym_var="EPSMEAN",bsym_var="bEPSMEAN"){
 temp_levels = c("No Symptoms","Symptoms")
 sym_var = paste(sym_var,cur_visit_id,sep=".")
 bsym_var = paste(bsym_var,cur_visit_id,sep=".")
 CATIE[,bsym_var] = factor(!is.na(CATIE[,sym_var])& CATIE[,sym_var]!=0)
 temp_unique = sort(unique(as.numeric(CATIE[,bsym_var])))
 levels(CATIE[,bsym_var]) = temp_levels[temp_unique ]
 CATIE[is.na(CATIE[,sym_var]),bsym_var] = NA
 CATIE[!is.na(CATIE[,sym_var])&CATIE[,sym_var]==0,sym_var] = NA
 return(CATIE)
}
###########################################################################
##  Several CATIE variables have upper and lower bounds, this is used for 
##   ensuring the bounds are met in the imputed data
###########################################################################
CATIE_get_bounds = function( CATIE, vars_list ){

 # set up the bounds
 upper_bounds = rep(NA,length(vars_list))
 names(upper_bounds) = vars_list
 lower_bounds = upper_bounds

 t = regexpr("BMI",vars_list)>0
 lower_bounds[t] = 15.0
 upper_bounds[t] = 65.0

 t_PANSS =  regexpr("PANSS",vars_list)
 upper_bounds[t_PANSS>0] = 210
 lower_bounds[t_PANSS>0] = 30

 # these are greater than 0
 # "MEDAD03","QOL_TOT", "YRS_PRES"
 for(jj in c("MEDAD03","QOL_TOT","YRS_PRES") ){
  t_temp =  regexpr(jj,vars_list,fixed=T)
  lower_bounds[ t_temp>0 ] = 0
 }
 t_temp =  regexpr("MEDAD03",vars_list)
 upper_bounds[t_temp>0] = 100
 
 t = regexpr("EPSMEAN",vars_list,fixed=T) > 0 & 
   regexpr("b",vars_list,fixed=T) < 0
 upper_bounds[t] = 4.0
 lower_bounds[t] = 0

 t = regexpr("QOL_TOT",vars_list,fixed=T) > 0
 upper_bounds[t] = 6.0

 return(list(upper_bounds=upper_bounds,lower_bounds=lower_bounds))
}

###############################################################################
####	Impute stage 2 treatment at the first follow-up visit	              
###############################################################################
CATIE_impute_stage2_treatment=function(CATIE){
 ## indicies of people in each of the randomization arms that have missing 
 ## stage 2 treatment
 temp_index_miss_stage2 = CATIE$SWITCH=="SWITCHED" & is.na(CATIE$TREATs2)
 temp_index_eff = temp_index_miss_stage2 & !is.na(CATIE$STAGE2_ARM) & 
 		  CATIE$STAGE2_ARM == "Efficacy"
 temp_index_tol = temp_index_miss_stage2 & !is.na(CATIE$STAGE2_ARM) & 
 		  CATIE$STAGE2_ARM == "Tolerability"

 Tol_t2 = c("Olanzapine","Risperidone","Ziprasidone")
 Eff_t2 = c("Olanzapine","Risperidone","Clozapine")
 ### if there is no one missing treatment then return CATIE as is
 if( sum(temp_index_miss_stage2)==0 ){
  print("No missing stage two treatment")
  return(CATIE)
 }
 ## for all possible treatment that a person could have been randomized 
 ## 	to before phase 2
 ## loop through so people cannot get randomized to the same treatment twice
 for( ii in  c("Olanzapine","Quetiapine","Risperidone","Perphenazine") ){

  ## These people recieved treatment ii in stage 1
  temp_treat_ii = ( CATIE[,"TREATs1"] == ii )
  
  if( ii != "Perphenazine" & ii != "Quetiapine" ){

   # for people who had either olanzapine or risperidone and diccontinued, make
   # sure not rerandomized to those medications again.
   ### Efficacy arm
   temp_treat_now = Eff_t2[-which(Eff_t2==ii)]   
   CATIE[temp_index_eff&temp_treat_ii,"TREATs2"] = temp_treat_now[
			   sample(size=sum(temp_treat_ii[temp_index_eff]),
			   x=1:2,replace=T,prob=c(.5,0.5)) ]
   ### Tolerability arm
   temp_treat_now = Tol_t2[-which(Tol_t2==ii)]   
   CATIE[temp_index_tol&temp_treat_ii,"TREATs2"] = temp_treat_now[
			   sample(size=sum(temp_treat_ii[temp_index_tol]),
			   x=1:2,replace=T,prob=c(.5,0.5)) ]
  }else{
   # As CATIE protocol randomization probabilities for clozapine 
   # and ziprasidone are higher than the other stage 2 treatments
   ### Efficacy arm
   temp_treat_now = Eff_t2
   CATIE[temp_index_eff&temp_treat_ii,"TREATs2"] = temp_treat_now[
			 sample(size=sum(temp_treat_ii[temp_index_eff]),
			 x=1:3,replace=T,prob=c(0.50,0.25,0.25)) ]
   ###              Tolerability arm
   temp_treat_now = Tol_t2
   CATIE[temp_index_tol&temp_treat_ii,"TREATs2"] = temp_treat_now[
			 sample(size=sum(temp_treat_ii[temp_index_tol]),
			 x=1:3,replace=T,prob=c(.50,0.25,0.25)) ]
   }
 }
 CATIE$TREAT.2 = CATIE$TREATs2
 CATIE$TREAT.2 = factor(CATIE$TREAT.2)
 levels(CATIE$TREAT.2) = c(levels(CATIE$TREAT.2),"Perphenazine","Quetiapine")
 CATIE$TREAT.2[CATIE$SWITCH=="STAYED"] = CATIE[CATIE$SWITCH=="STAYED",
 				       "TREATs1"]
 return(CATIE)
}

###########################################################################
## Make sure after imputations that binary indicator (bEPSMEAN) for no 
## symptoms matches with the continuous score (EPSMEAN)
###########################################################################
clean_up_side_effects_CATIE= function(CATIE, cur_visit_id, END){
 if(END){
  CATIE[!is.na(CATIE[,paste("END_bEPSMEAN",cur_visit_id,sep=".")])&
	  CATIE[,paste("END_bEPSMEAN",cur_visit_id,sep=".")]=="No Symptoms",
		   	paste("END_EPSMEAN",cur_visit_id,sep=".")] = 0
  CATIE[!is.na(CATIE[,paste("END_EPSMEAN",cur_visit_id,sep=".")])&
	  CATIE[,paste("END_EPSMEAN",cur_visit_id,sep=".")]==0,
	   	paste("END_bEPSMEAN",cur_visit_id,sep=".")] = "No Symptoms"
  }else{
   CATIE[!is.na(CATIE[,paste("bEPSMEAN",cur_visit_id,sep=".")])&
	CATIE[,paste("bEPSMEAN",cur_visit_id,sep=".")]=="No Symptoms",
			paste("EPSMEAN",cur_visit_id,sep=".")] = 0
   CATIE[!is.na(CATIE[,paste("EPSMEAN",cur_visit_id,sep=".")])&
	   CATIE[,paste("EPSMEAN",cur_visit_id,sep=".")]==0,
		paste("bEPSMEAN",cur_visit_id,sep=".")] = "No Symptoms"
  }
 return(CATIE)
}

##################################################################
### Impute missing PANSS variable with mixed effect model
##################################################################
CATIE_impute_PANSS=function( CATIE, cur_visit_id, id_var, END=TRUE,
			    time_varying_PANSS_predictors, SCH=TRUE,  
  	  	            variablesNoMiss, seed ){

 ## Keep track of how many visits have "gone by" 
 cur_visit_ids = 0:cur_visit_id

 # get all the time-varying variables with the visit id attached that will be
 #  used in the imputation. Either as predictors or "PANSSTOT" because it 
 #  needs to be imputed
 all_varying = NULL
 for(iv in cur_visit_ids){
  all_varying = c(all_varying,
  	      paste(c(time_varying_PANSS_predictors,"PANSSTOT"),iv,sep="." )) 
 }
 all_varying = unique(all_varying)
 variablesPredPANSSNoNest = setdiff(variablesNoMiss,all_varying)
 if( which(time_varying_PANSS_predictors =="MEDAD03") & 
     					 sum(names(CATIE)=="MEDAD03.0") == 0){
  CATIE$MEDAD03.0 = 0
 }
 if( which(time_varying_PANSS_predictors =="TREAT") & 
     					 sum(names(CATIE)=="TREAT.0") == 0){
  CATIE$TREAT.0 = "No Trt Yet"
  CATIE$TREAT.0 = factor(CATIE$TREAT.0)
 }
 # Separate out the data into that to be used in the imputation strategy 
 #    at this point and variables in the future
 CATIE_imp = CATIE[,unique(c(variablesPredPANSSNoNest,
		all_varying,"STAGE2_ARM"))]
 CATIE_other = CATIE[,unique( c(id_var, 
 	        setdiff(names(CATIE),names(CATIE_imp))))]
 ## reshape into long format to use the mixed effect imputation model.
 CATIE_long = reshape(CATIE_imp,direction="long", varying= all_varying,
			idvar=id_var,timevar="Visit")
 
 ### CATIE_long must be sorted by ID or pan will FREAK OUT
 CATIE_long = CATIE_long[sort(CATIE_long[,id_var],index.return=TRUE)$ix,]

 ## Create the design matrix to pass into pan
 temp_list = create_predictors4pan(CATIE_long,id_var=id_var,END=END,SCH=SCH,
			predictorsNoNest=variablesPredPANSSNoNest,
			predictorsNest=time_varying_PANSS_predictors)

 ##############################################################################
 ##############################################################################
 ####  This section is copied from comments in some example code from Schafers 
 ####  who created and maintains the pan function, please see his documentation
 ####  for more detailed information.
 ##############################################################################
 ## Sigma is the parameter that specifies the covariance matrix for response.
 ## Recall that the dimension of Sigma is (r x r) where r
 ## is the number of response variables (in this case, r=1). The prior
 ## distribution for Sigma is inverted Wishart with hyperparameters a 
 ## (scalar) and Binv (r x r), where a is the imaginary degrees of freedom
 ## and Binv/a is the prior guesstimate of Sigma. The value of a must be
 ## greater than or equal to r. The "least informative" prior possible
 ## would have a=r, so here we will take a=1. As a prior guesstimate of 
 ## Sigma we will use the (r x r) identity matrix, so Binv = 1*1 = 1.
 ##############################################################################
 ## By similar reasoning we choose the prior distribution for Psi. The
 ## dimension of Psi is (r*q x r*q) where q is the number of random
 ## effects in the model (which in this case is one).
 ## The hyperparameters for Psi are c and Dinv, where c is the
 ## imaginary degrees of freedom (which must be greater than or equal to
 ## r*q) and Cinv/d is the prior guesstimate of Psi. We will take d=1
 ## and Cinv=1*1 = 1.
 ##############################################################################
 y_pan = temp_list$resp
 subj_pan = as.integer(temp_list$subj)
 pred_pan = as.matrix(temp_list$Pred)
 attr(pred_pan,"contrasts") = NULL
 attr(pred_pan,"assign") = NULL
 xcol_pan = 1:ncol(pred_pan)
 zcol_pan = which(colnames(pred_pan)=="(Intercept)")
 prior_list=list(a=1, Binv=matrix(1,nrow=1,ncol=1), c=1, 
 		      Dinv=diag(1,nrow=1,ncol=1))
 pan_out = pan(y=y_pan,subj=subj_pan,pred=pred_pan, xcol=xcol_pan,
	       zcol=zcol_pan, prior=prior_list, seed=floor(runif(1)*32000), 
	       iter=250)
 # This is just to gaurd against collinearity, delete if fine at visit 2

 CATIE_long$PANSSTOT = pan_out$y  
 # make sure in bounds
 if( sum(CATIE_long$PANSSTOT<30) > 0 ){
  CATIE_long$PANSSTOT[CATIE_long$PANSSTOT<30] = 30
 }
 if( sum(CATIE_long$PANSSTOT>210) > 0 ){
  CATIE_long$PANSSTOT[CATIE_long$PANSSTOT<210] = 210
 }

 # Put the data back into wide format
 CATIE_impd = reshape(CATIE_long,direction="wide",idvar=id_var,timevar="Visit",
 	      v.names = c(time_varying_PANSS_predictors,"PANSSTOT"))
 CATIE_impd = CATIE_impd[sort(CATIE_impd[,id_var],index.return=TRUE)$ix,]
 CATIE_other = CATIE_other[sort(CATIE_other[,id_var],index.return=TRUE)$ix,]
 CATIE = merge(CATIE_other,CATIE_impd,by="CATIEID")

 if( sum(names(CATIE)=="MEDAD03.0")>0 ){
  CATIE = CATIE[,-which(names(CATIE)=="MEDAD03.0")]
 }
 if( sum(names(CATIE)=="TREAT.0")>0 ){
  CATIE = CATIE[,-which(names(CATIE)=="TREAT.0")]
 }
 return(CATIE)
}

###########################################################################
### This function makes everything numeric because pan needs all predictors
###   and outcomes numeric
###########################################################################
create_predictors4pan = function(CATIE_local,id_var,predictorsNoNest,
		      	              END=END,SCH=SCH,predictorsNest){
# This assumes there are no missing values in this data set other than PANSS
# ie that only the subjects that need to be imputed are in CATIE_local  

# predictorsNoNest = variablesPredPANSSNoNest
# predictorsNest = time_varying_PANSS_predictors

 ### with the results call pan such as
 # temp_list = create_predictors4pan(CATIE_local)
 # pan_out = pan(temp_list$resp,temp_list$subj,
 #             pred=temp_list$Pred, xcol=1:dim(X)[2] , zcol=1,
 #             prior=prior_list, seed=1,iter=1000);

 ### This is ordered categorical so just convert to numeric
 CATIE_pan = CATIE_local
 
 
 # if(sum(names(CATIE=="CS16.0"))>0){
 if( sum(regexpr(" ",names(CATIE_pan)) > 0)>0 ){
  CATIE_pan[,"CS16"] = as.numeric(CATIE_local[,"CS16"])  
 }
 ### Make Stage 2_Arm a variable we can use to nest end of stage predictors in 
 if(END & SCH){
  if( sum(regexpr("STAGE2_ARM",names(CATIE_pan)) > 0)>0 ){
   temp_levels = levels(CATIE_pan$STAGE2_ARM)
   CATIE_pan$STAGE2_ARM = as.numeric(CATIE_pan$STAGE2_ARM )
   CATIE_pan$STAGE2_ARM[is.na(CATIE_pan$STAGE2_ARM)] = 3
   CATIE_pan$STAGE2_ARM = factor(CATIE_pan$STAGE2_ARM,
   			  labels=c(temp_levels,"STAGE 1"))
  }
 }

 ### This is the imputation formula for PANSS
 ## predictorsNoNest - are those variabels not nested in end of phase
 ## predictorsNest - are those variables nested within SWITCH and 
 ## 		   STAGE2_ARM (which is not represented in one variable
 ##		   called STAGE2_ARM, just for this function).
 if(END){
  # if predicting end of stage PANSS nest predictors
  temp_form  = paste("~",paste(predictorsNoNest,collapse="+"),"+STAGE2_ARM*(",
  	          paste(predictorsNest,collapse="+"),")")
  CATIE_pan_mod = model.matrix(as.formula(temp_form), data=CATIE_pan)
 }else {
  # if only imputing scheduled PANSS no need to nest.
  temp_form  = paste("~",paste(predictorsNoNest,collapse="+")," + (",
  	          paste(predictorsNest,collapse="+"),")")
  CATIE_pan_mod = model.matrix(as.formula(temp_form), data=CATIE_pan)
 }

 ### return a list with things to input into pan
 ### Pred is matrix of predictors with no missing data, everythins is numeric
 ### 	   and it is in a matrix format
 ### response is PANSS in matrix form, with 1 column.
 ### subj is CATIEID as an integer.
 return(list(Pred=CATIE_pan_mod,
	resp=as.matrix(CATIE_local[,"PANSSTOT"],ncol=1),
	subj=as.integer(CATIE_local[,id_var])))
  
 ### with the results call pan such as
 # temp_list = create_predictors4pan(CATIE_local)
 # pan_out = pan(temp_list$resp,temp_list$subj,
 #             pred=temp_list$Pred, xcol=1:dim(X)[2] , zcol=1,
 #             prior=prior_list, seed=1, iter=1000);

}

###########################################################################
## Use an imputation model to impute if someone transitioned to the 
##     next treatment stage
###########################################################################
impute_CATIE_stage_transition = function( CATIE, variablesNoMiss, seed,
			      		  	 id_var=id_var ){

 ## Separate CATIE into data set involved in imputations,
 CATIE_imp = CATIE[,unique(c(variablesNoMiss,"SWITCH"))]
 CATIE_other = CATIE[,unique( c(id_var,
  	         setdiff(names(CATIE),c(variablesNoMiss,"SWITCH")) ))]

 ### Use MICE to impute missing scheduled variabels at this visit
 tmi = mice(data=CATIE_imp, m = 1, method="logreg", seed=seed, 
       			 printFlag=F, diagnostics=F, 
       			 defaultMethod=c("norm","logreg","polyreg") )
 CATIE_imp = complete(tmi)
 CATIE_imp = cbind(CATIE[,id_var],CATIE_imp)
 names(CATIE_imp)[1] = "CATIEID"
 CATIE = merge(CATIE_other, CATIE_imp)
 return(CATIE)

}

