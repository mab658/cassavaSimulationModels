if (!require("pacman")) {
  install.packages("pacman")
}
# Include all packages here
pacman::p_load(
  AlphaSimR,
  doParallel,
  foreach,
  RhpcBLASctl,
  rrBLUP,
  tidyverse,
  tidyr,
  tibble,
  here,
  prettycode,
  formattable
)
prettycode::prettycode()

setwd("/Users/mab658/Documents/cassavaSimulationModels")
here::set_here()


rm(list=ls())

# Load global parameters
source("./code/globalParameters.R")

# load a custom function to split population
# to split genetic resources for two mega-environments
source("./code/splitPop.R")

# load a function to generate simulated gxe data
# across multipe locations
source("./code/gxeSim.R")

# load a function to advance each mega environmente
# of narrow adaptation after burn-in
source("./code/megaEnv.R")

# Simulation runs for 10-year burn-in phase for nSimRun replicates
for(REP in 1:nSimRun){
  
  cat("\n Simulation run:",REP,"of", nSimRun,"replication", "\n")
  
  # Create founder genomes and initial parents population
  source("./code/createParents.R")
  
  # Fill breeding pipeline with unique indv. from initial parents
  source("./code/cassavaBreedingPipeline.R")
  
  # Run 10 years of burn-in by advancing years of breeding
  source("./code/burnin.R") # advance evaluation stages by years
}


# scenario 1 conventional breeding scheme
for(REP in 1:nSimRun){
  cat("\n Modeling scenario 1 conventional breeding scheme \n")

  cat("\n Simulation run:",REP,"of", nSimRun,"replication", "\n")

  # load the burnin phase for the second scenario
  load(paste0("./data/burnin_",REP,".rda"))

  # invoke the script to advance year of breeding for genomic selection
  source("./code/conv.R")
}

# Scenario 2: broad adaptation

for(REP in 1:nSimRun){
  
  cat("\n Modeling scenario 2 broad adapatation enabled GS \n")
  
  cat("\n Simulation run:",REP,"of", nSimRun,"replication", "\n")
  
  # load the burnin phase for the second scenario
  load(paste0("./data/burnin_",REP,".rda"))

  # invoke the script to advance year of breeding for genomic selection
  source("./code/broadAdaptation.R")
}


# Scenario 3: Narrow adaptation - megaEnv1
# This script implemented narrow  adaptation BP
# with GS to estimate BV in CET based on genomic data to boost 
# the selection accuracy at this early stage where phenotype 
# information is limited
# GS is used to advance from CET, PYT, AYT, and UYT

for(REP in 1:nSimRun){
  cat("\n Modeling scenario 3 narrow  adaptation enabled GS \n")
  
  cat("\n Simulation run:",REP,"of", nSimRun,"replication", "\n")

  # load the burnin phase for the second scenario
  load(paste0("./data/burnin_",REP,".rda"))
  
  # Train genomic selection (GS) model and safe the model fit
  gsModel <- RRBLUP(pop=trainPopME1, traits=1, use="pheno", simParam=SP)
  
  CET = setEBV(pop=CET, solution=gsModel)
  PYT = setEBV(pop=PYT, solution=gsModel)
  AYT = setEBV(pop=AYT, solution=gsModel)
  UYT = setEBV(pop=UYT, solution=gsModel)
  trainPop <- trainPopME1
  
  # Invoke megaEnv function to advance narrow-adaptation for megaEnv1
  output <- megaEnv(pval1UYT=0.6,pval2UYT=0.9,
                    pval1AYT=0.6,pval2AYT=0.7,
                    pvalPYT=0.7,pvalCET=0.5, 
                    pvalSDN=0.5)
  write.csv(output,file=paste0("./data/narrowME1_parms","_",REP,".csv"),row.names=FALSE)
}


# Scenario 3: Narrow adaptation - megaEnv2
# This script implemented narrow  adaptation BP
# with GS to estimate BV in CET based on genomic data to boost 
# the selection accuracy at this early stage where phenotype 
# information is limited
# GS is used to advance from CET, PYT, AYT, and UYT

for(REP in 1:nSimRun){
  cat("\n Modeling scenario 3 narrow  adaptation enabled GS \n")
  
  cat("\n Simulation run:",REP,"of", nSimRun,"replication", "\n")
  
  # load the burnin phase for the second scenario
  load(paste0("./data/burnin_",REP,".rda"))
  
  # Train genomic selection (GS) model and safe the model fit
  gsModel <- RRBLUP(pop=trainPopME2, traits=1, use="pheno", simParam=SP)
  
  # estimate gebv of the population
  CET = setEBV(pop=CET, solution=gsModel)
  PYT = setEBV(pop=PYT, solution=gsModel)
  AYT = setEBV(pop=AYT, solution=gsModel)
  UYT = setEBV(pop=UYT, solution=gsModel)
  
  trainPop <- trainPopME2
  output <- megaEnv(pval1UYT=0.1,pval2UYT=0.4,
                    pval1AYT=0.3,pval2AYT=0.4,
                    pvalPYT=0.3,pvalCET=0.5,
                    pvalSDN=0.5)
  write.csv(output,file=paste0("./data/narrowME2_parms","_",REP,".csv"),row.names=FALSE)
  #saveRDS(result,file = paste0("./data/narrowME2_parms","_",REP,".rds"))
}

(end.time <- Sys.time())
(end.time - start.time)