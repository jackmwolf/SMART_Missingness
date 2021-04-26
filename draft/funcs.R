source("MI_Shortreed.R")

############################################################################
#User inputs pct_missing, n_sims, n (later for future app) and mechanism

#################### Functions #######################
# Function to set up data for analysis after imputation
setup <- function(mi.dat) {
  # Delete P treatment
  d <- mi.dat
  mi.dat <- d[-which(d$TREAT.2 == "Perphenazine"), ]
  mi.dat <- mi.dat[-which(mi.dat$TREAT.1 == "Perphenazine"), ]
  n <- length(mi.dat$CATIEID)
  # Define dummy variables
  mi.dat$va <- ifelse(mi.dat$SITETYPE == "VA", 1, 0)
  mi.dat$combo <- ifelse(mi.dat$SITETYPE == "Combo", 1, 0)
  mi.dat$uni <- ifelse(mi.dat$SITETYPE == "University Clinic", 1, 0)
  mi.dat$state <- ifelse(mi.dat$SITETYPE == "State Mental Health", 1, 0)
  mi.dat$EMPLOY <- ifelse(mi.dat$EMPLOY == "Did Not Work", 0, 1)
  mi.dat$EXACER <- ifelse(mi.dat$EXACER == "Yes", 1, 0)
  mi.dat$SEX <- ifelse(mi.dat$SEX == "Male", 1, 0)
  mi.dat$TD <- ifelse(mi.dat$TD == "Yes", 1, 0)
  # Treatment 1 dummies
  mi.dat$qtrt1 <- ifelse(mi.dat$TREAT.1 == "Quetiapine", 1, 0)
  mi.dat$otrt1 <- ifelse(mi.dat$TREAT.1 == "Olanzapine", 1, 0)
  mi.dat$rtrt1 <- ifelse(mi.dat$TREAT.1 == "Risperidone", 1, 0)
  mi.dat$ptrt1 <- ifelse(mi.dat$TREAT.1 == "Perphenazine", 1, 0)
  # Treatment 2 dummies
  mi.dat$qtrt2 <- ifelse(mi.dat$TREAT.2 == "Quetiapine", 1, 0)
  mi.dat$otrt2 <- ifelse(mi.dat$TREAT.2 == "Olanzapine", 1, 0)
  mi.dat$ctrt2 <- ifelse(mi.dat$TREAT.2 == "Clozapine", 1, 0)
  mi.dat$rtrt2 <- ifelse(mi.dat$TREAT.1 == "Risperidone", 1, 0)
  mi.dat$ztrt2 <- ifelse(mi.dat$TREAT.2 == "Ziprasidone", 1, 0)
  mi.dat$SWITCH <- ifelse(mi.dat$SWITCH == "SWITCHED", 1, 0)
  mi.dat$STAGE2_ARM <- ifelse(mi.dat$STAGE2_ARM == "Efficacy", 1, 0) # 1 is efficacy
  mi.dat$pan <- (mi.dat$PANSSTOT.0 + mi.dat$PANSSTOT.1 + mi.dat$PANSSTOT.2) / 3
  levels(mi.dat$TREAT.2)[levels(mi.dat$TREAT.2) == "Ziprasidone" | levels(mi.dat$TREAT.2) == "Clozapine"] <- "Clz"
  return(mi.dat)
}

#Function to induce missingness
gen.mis.mnar = function(N,dropout.time,psi.s1,psi.s2,psi.s3,dat){
  set.seed(1996)
  r1pred = exp(psi.s1[1]+psi.s1[2]*dat$BMI.1+psi.s1[3]*dat$SWITCH+
                 psi.s1[4]*dat$MEDAD03.1+psi.s1[5]*dat$PANSSTOT.1) 
  r2pred = exp(psi.s2[1]+psi.s2[2]*dat$SWITCH+
                 psi.s2[3]*dat$MEDAD03.1+psi.s2[4]*dat$PANSSTOT.1) 
  r3pred = exp(psi.s3[1]+psi.s3[2]*dat$BMI.2+
                 psi.s3[3]*dat$MEDAD03.2+psi.s3[4]*dat$PANSSTOT.2) 
  #No treatment assignment
  pr.r1 <-  r1pred/(r1pred+1) 
  if(sum(is.na(pr.r1))>0){pr.r1[which(is.na(pr.r1))]=1};summary(pr.r1)
  mis.s1= rbinom(N,1,pr.r1);sum(mis.s1==1)
  
  #Stage 1 + no treatment assignment
  pr.r2 <-  r2pred/(r2pred+1); summary(pr.r2)
  mis.s2= rbinom(N,1,pr.r2);sum(mis.s2==1)
  
  #Stage 2
  pr.r3 <-  r3pred/(r3pred+1);summary(pr.r3)
  mis.s3= rbinom(N,1,pr.r3);sum(mis.s3==1)
  
  #Induce missingness
  dat.mi=dat;nvar=length(names(dat.mi))
  switch.i = which(names(dat.mi)=="SWITCH")
  bmi1.i = which(names(dat.mi)=="BMI.1")
  bmi2.i = which(names(dat.mi)=="BMI.2")
  
  if(dropout.time==1){
    dat.mi$YRS_PRES[runif(8,1,length(dat.mi$CATIEID))]=NA
    dat.mi[which(mis.s1==1),c(switch.i:nvar)]=NA
    dat.mi[which(mis.s2==1),c(switch.i:(bmi1.i-1))]=NA
    dat.mi[which(mis.s2==1),c(bmi2.i:nvar)]=NA
    dat.mi[which(mis.s3==1),c(bmi2.i:nvar)]=NA}  
  if(dropout.time==2){
    dat.mi$YRS_PRES[runif(8,1,length(dat.mi$CATIEID))]=NA
    dat.mi[which(mis.s2==1),c(switch.i:bmi1.i-1)]=NA
    dat.mi[which(mis.s2==1),c(bmi2.i:nvar)]=NA
    dat.mi[which(mis.s3==1),c(bmi2.i:nvar)]=NA}
  if(dropout.time==3){
    dat.mi$YRS_PRES[runif(8,1,length(dat.mi$CATIEID))]=NA
    dat.mi[which(mis.s3==1),c(bmi2.i:nvar)]=NA}
  dat.mi$SWITCH[which(dat.mi$SWITCH==1)]="SWITCHED"
  dat.mi$SWITCH[which(dat.mi$SWITCH==0)]="STAYED"
  dat.mi$SWITCH=as.factor(dat.mi$SWITCH)
  dat.mi$bEPSMEAN.0=NA
  return(dat.mi)}

#Function to generate MAR 
gen.mis.mar = function(N,dropout.time,psi.s1,psi.s2,psi.s3,dat){
  set.seed(1996)
  r1pred = exp(psi.s1[1]+psi.s1[2]*dat$BMI.0+psi.s1[3]*dat$PANSSTOT.0) 
  r2pred = exp(psi.s2[1]+psi.s2[2]*dat$MEDAD03.1+psi.s2[3]*dat$PANSSTOT.1) 
  r3pred = exp(psi.s3[1]+psi.s3[2]*dat$BMI.1+ psi.s3[3]*dat$SWITCH+
                 psi.s3[4]*dat$MEDAD03.1+psi.s3[5]*dat$PANSSTOT.1) 
  #No treatment assignment
  pr.r1 <-  r1pred/(r1pred+1) 
  if(sum(is.na(pr.r1))>0){pr.r1[which(is.na(pr.r1))]=1};summary(pr.r1)
  mis.s1= rbinom(N,1,pr.r1);sum(mis.s1==1)
  
  #Stage 1 + no treatment assignment
  pr.r2 <-  r2pred/(r2pred+1); summary(pr.r2)
  mis.s2= rbinom(N,1,pr.r2);sum(mis.s2==1)
  
  #Stage 2
  pr.r3 <-  r3pred/(r3pred+1);summary(pr.r3)
  mis.s3= rbinom(N,1,pr.r3);sum(mis.s3==1)
  
  #Induce missingness
  dat.mi=dat;nvar=length(names(dat.mi))
  switch.i = which(names(dat.mi)=="SWITCH")
  bmi1.i = which(names(dat.mi)=="BMI.1")
  bmi2.i = which(names(dat.mi)=="BMI.2")
  
  if(dropout.time==1){
    dat.mi$YRS_PRES[runif(8,1,length(dat.mi$CATIEID))]=NA
    dat.mi[which(mis.s1==1),c(switch.i:nvar)]=NA
    dat.mi[which(mis.s2==1),c(switch.i:(bmi1.i-1))]=NA
    dat.mi[which(mis.s2==1),c(bmi2.i:nvar)]=NA
    dat.mi[which(mis.s3==1),c(bmi2.i:nvar)]=NA}  
  if(dropout.time==2){
    dat.mi$YRS_PRES[runif(8,1,length(dat.mi$CATIEID))]=NA
    dat.mi[which(mis.s2==1),c(switch.i:bmi1.i-1)]=NA
    dat.mi[which(mis.s2==1),c(bmi2.i:nvar)]=NA
    dat.mi[which(mis.s3==1),c(bmi2.i:nvar)]=NA}
  if(dropout.time==3){
    dat.mi$YRS_PRES[runif(8,1,length(dat.mi$CATIEID))]=NA
    dat.mi[which(mis.s3==1),c(bmi2.i:nvar)]=NA}
  dat.mi$SWITCH[which(dat.mi$SWITCH==1)]="SWITCHED"
  dat.mi$SWITCH[which(dat.mi$SWITCH==0)]="STAYED"
  dat.mi$SWITCH=as.factor(dat.mi$SWITCH)
  dat.mi$bEPSMEAN.0=NA
  return(dat.mi)}

#Function to generate MCAR missingness
gen.mis.mcar = function(N,dropout.time,pr.r1, pr.r2, pr.r3,dat){
  set.seed(1996)
  #No treatment assignment
  pr.r1 <-  .25
  if(sum(is.na(pr.r1))>0){pr.r1[which(is.na(pr.r1))]=1};summary(pr.r1)
  mis.s1= rbinom(N,1,pr.r1);sum(mis.s1==1)
  #Stage 1 + no treatment assignment
  pr.r2 <-  .1; summary(pr.r2)
  mis.s2= rbinom(N,1,pr.r2);sum(mis.s2==1)
  #Stage 2
  pr.r3 <-  .25;summary(pr.r3)
  mis.s3= rbinom(N,1,pr.r3);sum(mis.s3==1)
  
  #Induce missingness
  dat.mi=dat;nvar=length(names(dat.mi))
  switch.i = which(names(dat.mi)=="SWITCH")
  bmi1.i = which(names(dat.mi)=="BMI.1")
  bmi2.i = which(names(dat.mi)=="BMI.2")
  
  if(dropout.time==1){
    dat.mi$YRS_PRES[runif(8,1,length(dat.mi$CATIEID))]=NA
    dat.mi[which(mis.s1==1),c(switch.i:nvar)]=NA
    dat.mi[which(mis.s2==1),c(switch.i:(bmi1.i-1))]=NA
    dat.mi[which(mis.s2==1),c(bmi2.i:nvar)]=NA
    dat.mi[which(mis.s3==1),c(bmi2.i:nvar)]=NA}  
  if(dropout.time==2){
    dat.mi$YRS_PRES[runif(8,1,length(dat.mi$CATIEID))]=NA
    dat.mi[which(mis.s2==1),c(switch.i:bmi1.i-1)]=NA
    dat.mi[which(mis.s2==1),c(bmi2.i:nvar)]=NA
    dat.mi[which(mis.s3==1),c(bmi2.i:nvar)]=NA}
  if(dropout.time==3){
    dat.mi$YRS_PRES[runif(8,1,length(dat.mi$CATIEID))]=NA
    dat.mi[which(mis.s3==1),c(bmi2.i:nvar)]=NA}
  dat.mi$SWITCH[which(dat.mi$SWITCH==1)]="SWITCHED"
  dat.mi$SWITCH[which(dat.mi$SWITCH==0)]="STAYED"
  dat.mi$SWITCH=as.factor(dat.mi$SWITCH)
  dat.mi$bEPSMEAN.0=NA
  return(dat.mi)}


regime.mean=function(trt1,trt2,dat){
  p=.5 #Randomization probability 
  if(trt1=="Quetiapine"){
    if(trt2=="Olanzapine"){p=.25}
    if(trt2=="Risperidone"){p=.25}
    if(trt2=="Clz"){p=.5}}
  reg1w=ifelse(dat$TREAT.1==trt1,1/.33,0)
  reg1w=ifelse(dat$TREAT.1==trt1 & dat$TREAT.2==trt2,reg1w*(1/p),reg1w)
  weighted.mean(dat$pan,reg1w)}

#Function for running the bootstrap 
regime.mean.boot=function(data,indices,trt1,trt2){
  data.boot <- data[indices,]
  p=.5
  if(trt1=="Quetiapine"){
    if(trt2=="Olanzapine"){p=.25}
    if(trt2=="Risperidone"){p=.25}
    if(trt2=="Clz"){p=.5}}
  reg1w=ifelse(data.boot$TREAT.1==trt1,1/.33,0)
  reg1w=ifelse(data.boot$TREAT.1==trt1 & data.boot$TREAT.2==trt2,reg1w*(1/p),reg1w)
  weighted.mean(data.boot$pan,reg1w)}

#Function to find SE from multiple imputed sets (Rubin's rule)
se.mi = function(reg, var){
  B = var(reg)
  W = mean(var)
  T = W + (1+(1/5))*B
  return(sqrt(T))}

#True means of population data 
true.means <- c(75.79, 75.73, 76.03, 76.06, 78.50, 77.90, 77.23)
t.regO1 <- 75.79; t.regO2 <- 75.73; t.regR1 <- 76.03; t.regR2 <- 76.06 
t.regQ1 <- 78.50; t.regQ2 <- 77.90; t.regQ3 <- 77.23

#Main function
main <- function(sims, pct_mis, mis_mech, samp_size, true.dat) {
  set.seed(sims)
    # Take a sample of size n from the population
    sim.dat <- sample_n(true.dat, samp_size, replace = FALSE)
    
    #### Induce missingness depending on pct_missing and missingness mechanism ####
    if (mis_mech=="MNAR" & pct_mis=="60"){
      psi.s1 <- c(100,6,30,-10,5) #BMI,SWITCH,MEDAD,PANSS
      psi.s2 <- c(100,80,-10,5) #SWITCH,MEDAD,PANSS
      psi.s3 <- c(20,2,-10,6) #BMI,MEDAD,PANSS 
      #Induce missingness in dataset
      dat.mis <- gen.mis.mnar(N=samp_size,dropout.time=1,psi.s1=psi.s1,psi.s2 = psi.s2,psi.s3=psi.s3,dat=sim.dat)
    }
    
    if (mis_mech=="MNAR" & pct_mis=="30"){
      psi.s1 <- c(30,6,30,-20,5) #BMI,SWITCH,MEDAD,PANSS
      psi.s2 <- c(42,80,-15,5) #SWITCH,MEDAD,PANSS
      psi.s3 <- c(45,2,-10,6) #BMI,MEDAD,PANSS 
      #Induce missingness
      dat.mis <- gen.mis.mnar(N=samp_size,dropout.time=1,psi.s1=psi.s1,psi.s2 = psi.s2,psi.s3=psi.s3,dat=sim.dat)
    }
  
    if (mis_mech=="MAR" & pct_mis==60){
      psi.s1 <- c(-6.2,.05,.05) #Intercept, baseline BMI, PANSS
      psi.s2 <- c(100,-10,5) #Intercept, baseline MEDAD, PANSS.1
      psi.s3 <- c(20,2,80,-10,6) #Intercept, BMI.1, SWITCH, MEDAD.1, PANSS.1
      #Induce missingness
      dat.mis <- gen.mis.mar(N=samp_size,dropout.time=1,psi.s1=psi.s1,psi.s2 = psi.s2,psi.s3=psi.s3,dat=sim.dat)
    }
    
    if (mis_mech=="MAR" & pct_mis==30){
      psi.s1 <- c(-8,.05,.05) #Intercept, baseline BMI, PANSS
      psi.s2 <- c(100,-10,5) #Intercept, baseline MEDAD, PANSS.1
      psi.s3 <- c(20,2,80,-11,6) #Intercept, BMI.1, SWITCH, MEDAD.1, PANSS.1
      #Induce missingness
      dat.mis <- gen.mis.mar(N=samp_size,dropout.time=1,psi.s1=psi.s1,psi.s2 = psi.s2,psi.s3=psi.s3,dat=sim.dat)
    }
    
    if (mis_mech =="MCAR" & pct_mis == 60){
      dat.mis <- gen.mis.mcar(N=samp_size,dropout.time=1,pr.r1=.25, pr.r2=.1, pr.r3=.25,dat=sim.dat)
    }
    
    if (mis_mech =="MCAR" & pct_mis == 30){
      dat.mis <- gen.mis.mcar(N=samp_size,dropout.time=1,pr.r1=.15, pr.r2=.05, pr.r3=.15,dat=sim.dat)
    }
    
    #### Perform MI ####
    varO1 = varO2 = varR1 = varR2 = varQ1 = varQ2 = varQ3 <- c()
    regO1 = regO2 = regR1 = regR2 = regQ1 = regQ2 = regQ3 <- c()
    sim.O1 = sim.O2 = sim.R1 = sim.R2 = sim.Q1 = sim.Q2 = sim.Q3 = 0
    CovO1 = CovO2 = CovR1 = CovR2 = CovQ1 = CovQ2 = CovQ3 = 0
    #Run MI 10 times and get 10 estimates
    for (i in 1:10){
      set.seed(i)
      imp.dat = smart.mi(dat=dat.mis,mis.bas=FALSE,mis.end=TRUE,mis.switch=TRUE)
      
      #Set up data
      ipw.dat = setup(mi.dat=imp.dat)
      
    #### Regime means #### 
      #Regime means and SE's via IPW
      # R is the # of bootstrap simulations
      # (oo,or), (oo,cl or z), (qq,qo), (qq,qr), (qq,qcl or qz), (rr,ro), (rr,rcl or rz)#
      regO1[i] = regime.mean(trt1="Olanzapine",trt2="Risperidone",dat=ipw.dat)
      reg1.boot = boot(data=ipw.dat,statistic=regime.mean.boot,trt1="Olanzapine",trt2="Risperidone",R=1000)
      varO1[i] = var(reg1.boot$t)
      
      regO2[i] = regime.mean(trt1="Olanzapine",trt2="Clz",dat=ipw.dat)
      reg2.boot = boot(data=ipw.dat,statistic=regime.mean.boot,trt1="Olanzapine",trt2="Clz",R=1000)
      varO2[i] = var(reg2.boot$t)
      
      regR1[i] = regime.mean(trt1="Risperidone",trt2="Olanzapine",dat=ipw.dat)
      reg3.boot = boot(data=ipw.dat,statistic=regime.mean.boot,trt1="Risperidone",trt2="Olanzapine",R=1000)
      varR1[i] = var(reg3.boot$t)
      
      regR2[i] = regime.mean(trt1="Risperidone",trt2="Clz",dat=ipw.dat)
      reg4.boot = boot(data=ipw.dat,statistic=regime.mean.boot,trt1="Risperidone",trt2="Clz",R=1000)
      varR2[i] = var(reg4.boot$t)
      
      regQ1[i] = regime.mean(trt1="Quetiapine",trt2="Olanzapine",dat=ipw.dat)
      reg5.boot = boot(data=ipw.dat,statistic=regime.mean.boot,trt1="Quetiapine",trt2="Olanzapine",R=1000)
      varQ1[i] = var(reg5.boot$t)
      
      regQ2[i] = regime.mean(trt1="Quetiapine",trt2="Risperidone",dat=ipw.dat)
      reg6.boot = boot(data=ipw.dat,statistic=regime.mean.boot,trt1="Quetiapine",trt2="Risperidone",R=1000)
      varQ2[i] = var(reg6.boot$t)
      
      regQ3[i]= regime.mean(trt1="Quetiapine",trt2="Clz",dat=ipw.dat)
      reg7.boot = boot(data=ipw.dat,statistic=regime.mean.boot,trt1="Quetiapine",trt2="Clz",R=1000)
      varQ3[i] = var(reg7.boot$t) 
    }
    
    #### Final results from simulation ####
    #Regime Means
    sim.O1 <- mean(regO1); sim.O2 <- mean(regO2); sim.R1 <- mean(regR1); sim.R2 <- mean(regR2)
    sim.Q1 <- mean(regQ1); sim.Q2 <- mean(regQ2); sim.Q3 <- mean(regQ3)
    
    #Regime SE's
    sim.seO1 <- se.mi(reg=regO1,var=varO1); sim.seO2 <- se.mi(reg=regO2,var=varO2); 
    sim.seR1 <- se.mi(reg=regR1,var=varR1); sim.seR2 <- se.mi(reg=regR2,var=varR2)
    sim.seQ1 <- se.mi(reg=regQ1,var=varQ1); sim.seQ2 <- se.mi(reg=regQ2,var=varQ2); 
    sim.seQ3 <- se.mi(reg=regQ3,var=varQ3)
    
    #CI and Coverage probability
    CIO1 <- sim.O1+c(-1,1)*qnorm(.975)*sim.seO1; CovO1 <- ifelse(t.regO1>CIO1[1] & t.regO1<CIO1[2],1,0)
    CIO2 <- sim.O2+c(-1,1)*qnorm(.975)*sim.seO2; CovO2 <- ifelse(t.regO2>CIO2[1] & t.regO2<CIO2[2],1,0)
    CIR1 <- sim.R1+c(-1,1)*qnorm(.975)*sim.seR1; CovR1 <- ifelse(t.regR1>CIR1[1] & t.regR1<CIR1[2],1,0)
    CIR2 <- sim.R2+c(-1,1)*qnorm(.975)*sim.seR2; CovR2 <- ifelse(t.regR2>CIR2[1] & t.regR2<CIR2[2],1,0)
    CIQ1 <- sim.Q1+c(-1,1)*qnorm(.975)*sim.seQ1; CovQ1 <- ifelse(t.regQ1>CIQ1[1] & t.regQ1<CIQ1[2],1,0)
    CIQ2 <- sim.Q2+c(-1,1)*qnorm(.975)*sim.seQ2; CovQ2 <- ifelse(t.regQ2>CIQ2[1] & t.regQ2<CIQ2[2],1,0)
    CIQ3 <- sim.Q3+c(-1,1)*qnorm(.975)*sim.seQ3; CovQ3 <- ifelse(t.regQ3>CIQ3[1] & t.regQ3<CIQ3[2],1,0)
    
    sim.means <- c(sim.O1,sim.O2,sim.R1,sim.R2,sim.Q1,sim.Q2,sim.Q3)
    sim.cov <- c(CovO1,CovO2,CovR1,CovR2,CovQ1,CovQ2,CovQ3)
    return(cbind(sim.means,sim.cov))
  }
################################################################################
# 
# #Aggregate over n_sims
# sim.res <- data.frame(t(sapply(1:n_sims, main, pct_mis = pct_missing, mis_mec = mechanism, samp_size = 1000)))
# 
# ########################## Final Results #######################################
# 
# res <- sapply(sim.res, mean)
# names(res)=c("O1","O2","R1","R2","Q1","Q2","Q3","C.O1","C.O2","C.R1","C.R2","C.Q1","C.Q2","C.Q3")
# 
# #Metrics send to user
# bias <- res[1:7]-true.means
# se <- sapply(sim.res[1:7],sd)
# cov.probability <- res[8:14]
# mse <- c()
# for (i in 1:7){
#   m <- mean( (sim.res[,i] - true.means[i])^2 )
#   mse <- c(m,mse)}
# mse <- round(mse)
# 
# 
# #Plot
# Regime <- c(rep("DTR1",n_sims),rep("DTR2",n_sims),rep("DTR3",n_sims),
#             rep("DTR4",n_sims),rep("DTR5",n_sims),rep("DTR6",n_sims),
#             rep("DTR7",n_sims))
# 
# reg.means <- c(sim.res[,1],sim.res[,2],sim.res[,3],sim.res[,4],sim.res[,5],
#                sim.res[,6],sim.res[,7])
# 
# plot.df <- data.frame(Regime, reg.means)
# 
# plot = ggplot(plot.df, aes(Regime, reg.means, fill=Regime)) +
#   geom_boxplot() +
#   labs(y = "Expected PANSS score", x="Embedded DTR",
#        title="Comparing embedded DTRs")
# 
# #Add the true means to the plot (red dots)
# plot +
#   annotate("point", x = "DTR1", y = true.means[1], colour = "red",size=3)+
#   annotate("point", x = "DTR2", y = true.means[2], colour = "red",size=3)+
#   annotate("point", x = "DTR3", y = true.means[3], colour = "red",size=3)+
#   annotate("point", x = "DTR4", y = true.means[4], colour = "red",size=3)+
#   annotate("point", x = "DTR5", y = true.means[5], colour = "red",size=3)+
#   annotate("point", x = "DTR6", y = true.means[6], colour = "red",size=3)+
#   annotate("point", x = "DTR7", y = true.means[7], colour = "red",size=3)+
#   theme(legend.position = "none")
# 
