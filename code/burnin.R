# After populating the cassava breeding pipeline
# we implement the burn-in phase by advancing population of individuals
# through the breeding scheme over years
# Burnin phase serves as common starting point for
# different breeding scenarios (programs)
# Advance breeding by year by working backwards through
# pipeline  to avoid copying data

# we track the progress of breeding program
# by initializing variables to save the genetic parameters
simRun <- vector("numeric",length=nCycles)
meanGV <- vector("numeric",length=nCycles)
varGV <- vector("numeric",length=nCycles)

accPheno <- vector("numeric",length=nCycles)
accEbv <- vector("numeric",length=nCycles)
result <- data.frame()

#prentMean <- parentVar <- PYTMean <- PYTVar <-  matrix(NA, nrow = nCycles, ncol = 1)

# p-values for breeding cycles in burnin and future breeding phase
#P = runif(nCycles)

# Run 10 years of burn-in (Cycle years) as a common starting point
#- for different scenarios.
# This marks the beginning of a new breeding cycle where new parents are
# formed and taking to crossing block

# Advance breeding program by Working backwards through pipeline
# to avoid copying data

for(year in 1:burninYears){#Change to any number of desired years
  cat("Burnin phase cycle year:",year,"of",burninYears, "\n")

  # Select variety to be released

  variety <- selectInd(pop=UYT,nInd=nVarietySel, use="pheno")

  # Year 5 Uniform Yield Trial (UYT)

  UYT <- selectInd(pop=AYT, nInd=nUYT, use="pheno")
  UYT <- setPheno(pop=UYT, varE=errVarUYT, reps=repUYT, simParam=SP)

  # Year 4 Advance Yield Trial (AYT)
  AYT <- selectInd(pop=PYT, nInd = nAYT, use = "pheno")
  AYT <- setPheno(pop=AYT, varE=errVarAYT, reps=repAYT,simParam=SP)

  # Year 3 Preliminary  Yield  Trial (PYT)
  PYT <- selectInd(pop=CET, nInd=nPYT, use="pheno", simParam=SP)
  PYT <- setPheno(pop=PYT, varE=errVarPYT, reps=repPYT,simParam=SP)


  # Year 2 Clonal Evaluatiion Trial (CET)
  CET <- selectWithinFam(pop=SDN, nInd=famSize, use="pheno")
  CET <- selectInd(pop=CET, nInd=nCET, use="pheno", simParam=SP)
  CET <-  setPheno(pop=CET, varE=errVarCET, reps=repCET,simParam=SP)

  # Year 1  Seedling Nursery
  SDN <- F1
  SDN <- setPheno(pop=SDN, varE=errVarSDN,reps=repSDN,simParam=SP)

  # recycle new parents in each year of burn-in phase
  # by selecting individuals from joint (AYT, PYT & CET) population
  parents <- selectInd(c(UYT,AYT,PYT), nInd=50, use="pheno")

  # Year 1  crossing
  F1 <- randCross(pop=parents, nCrosses=nCrosses,
                  nProgeny=nProgeny, simParam=SP)


  # Constitute the training population (TP)
  # from 8-year cycle to 10-year burn-in phase
  # using CET,PYT, and AYT population to have 3 years of TP population genotyped

  # NOTE: Training records are collected from CET,PYT,and AYT with
  # a sliding-window process.
  # Accumulation of the records starts in the burn-in year 'startTrainPop'.
  # Once the burn-in period is over, the sliding-window process removes the oldest
  # records. The number of years retained is set with "limityear" in GlobalParameters.

  if(year == startTrainPop){
    trainPop = c(CET, PYT, AYT)
  } else if(year > startTrainPop && year <= burninYears){
    trainPop = c(trainPop, CET, PYT, AYT)
  }

  meanGV[year] <- meanG(PYT)
  varGV[year] <- varG(PYT)
  }

gGain <- meanGV - meanGV[burninYears]
relVar <- varGV/varGV[burninYears]

# Save burn-in phase as global environment for each simulation run
# to  be loaded later for other scenarios
# Save the BURNIN
cat(paste0("saving burn-in for REP ", REP, "\n"))
save.image(paste0("../data/burnin_",REP,".rda"))
