# After populating the cassava breeding pipeline
# we implement the burn-in phase by advancing population of individuals
# through the breeding scheme over years where the parents are recycled.
# Burnin phase serves as common starting point for
# different breeding scenarios (programs)
# Advance breeding by year by working backwards through
# pipeline  to avoid copying data

# we track the progress of breeding program
# by initializing variables to save the genetic parameters
simRun <- vector("numeric",length=nCycles)
meanGV <- vector("numeric",length=nCycles)
varGV <- vector("numeric",length=nCycles)

accCET <- accPYT <- accAYT <- accUYT <- vector("numeric",length=nCycles)

# Run 10 years of burn-in (Cycle years) as a common starting point
# for different scenarios.
# This marks the beginning of a new breeding cycle where new parents are
# formed and taking to crossing block

# Advance breeding program by Working backwards through pipeline
# to avoid copying data

for(year in 1:burninYears){#Change to any number of desired years
  cat("Burnin phase cycle year:",year,"of",burninYears, "\n")

  # Select variety to be released

  variety <- selectInd(pop=UYT,nInd=nVarietySel, use="pheno",simParam=SP)

  # Uniform Yield Trial (UYT)

  UYT <- selectInd(pop=AYT, nInd=nUYT, use="pheno",simParam=SP)
  UYT <- gxeSim(pval1=0.1, pval2=0.9, pop=UYT,
                varE=errVarUYT, nreps=repUYT,
                locSize=8)
  

  # Advance Yield Trial (AYT)
  AYT <- selectInd(pop=PYT, nInd=nAYT, use="pheno",simParam=SP)
  AYT <- gxeSim(pval1=0.3, pval2=0.7, pop=AYT,
                varE=errVarAYT, nreps=repAYT,locSize=4)
  

  # Preliminary  Yield  Trial (PYT)
  PYT <- selectInd(pop=CET, nInd=nPYT, use="pheno", simParam=SP)
  PYT <- gxeSim(pval1=0.3, pval2=0.7, pop=PYT,
                varE=errVarPYT,nreps=repPYT,locSize=2)


  #  Clonal Evaluatiion Trial (CET)
  CET <- selectWithinFam(pop=SDN, nInd=famSize, use="pheno", simParam=SP)
  CET <- selectInd(pop=CET, nInd=nCET, use="pheno", simParam=SP)
  CET <-  setPheno(pop=CET, varE=errVarCET, reps=repCET,p=0.5,simParam=SP)

  # Seedling Nursery
  SDN <- setPheno(pop=F1, varE=errVarSDN,reps=repSDN,
                  p=0.5, simParam=SP)

  # recycle new parents in each year of burn-in phase
  # by selecting best individuals from joint (UYT and AYT) population
  parents <- selectInd(c(UYT,AYT), nInd=50, 
                       use="pheno",simParam=SP)

  #  crossing block
  F1 <- randCross(pop=parents, nCrosses=nCrosses,
                  nProgeny=nProgeny, simParam=SP)
  
  # report the mean and variance of estimated GV in each year
  meanGV[year] <- meanG(PYT)
  varGV[year] <- varG(PYT)
  
  # Constitute the training population (TP)
  # from 8-year cycle to 10-year of burn-in phase
  # using CET,PYT, AYT, and UYT population to have 3 years of TP

  # NOTE: Training records are collected from CET,PYT,AYT,and UYT with
  # a sliding-window process.
  # Accumulation of the records starts in the burn-in year 'startTrainPop'.
  # Once the burn-in period is over, the sliding-window process removes the oldest
  # records. The number of years retained is set with "limityear" in GlobalParameters.

  if(year == startTrainPop){
    trainPop = c(CET, PYT, AYT, UYT) # TP for broad adaptation

    # Split the CET to constitute TP for mega-environments- ME1 and ME2 in narrow-adaptation
    splitCET <- splitPop(CET)

    trainPopME1 <- c(splitCET$ME1,PYT, AYT, UYT)
    trainPopME2 <- c(splitCET$ME2,PYT, AYT, UYT)
  } else if(year > startTrainPop && year <= burninYears){
    
    # update TP from year 9 to  10 of burnin phase
    trainPop = c(trainPop, CET, PYT, AYT, UYT) # TP broad adaptation

    splitCET <- splitPop(CET) # split TP set for ME1 and ME2

    trainPopME1 <- c(trainPopME1,splitCET$ME1,PYT, AYT, UYT)
    trainPopME2 <- c(trainPopME2,splitCET$ME2,PYT, AYT, UYT)
  }
} # end 10-year burn-in

# Save the state of simulation at final year=10 of  burn-in phase as global environment 
# for each simulation replicate. This will be loaded for other scenarios

cat(paste0("saving burn-in for REP ", REP, "\n"))
save.image(paste0("./data/burnin_",REP,".rda"))
