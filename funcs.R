source("MI_Shortreed.R")


main <- function(n, pct_missing, mechanism, n_sims) {
  
  
  # Aggregate over n_sims repititions
  {
    # Take a sample of size n from the population
    sim.dat <- sample_n(true.dat, n, replace = FALSE)
    
    # Induce missingness
    
    # Perform MI
    
  }
  
}



# Function to set up data
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
