# Scenario 1
# This script implemented broad adaptation  breeding programs
# with GS to estimate BV in CET based on genomic data to boost 
# the selection accuracy at this early stage where phenotype 
# information is limited
# GS is used to advance from CET, PYT, AYT, and UYT


# Train genomic selection (GS) model and safe the model fit
gsModel <- RRBLUP(pop=trainPop, traits=1, use="pheno", simParam=SP)

# Estimate BVs for CET, PYT, AYT, and UYT for being used to
# constitute the training population set and they are genotyped
CET = setEBV(pop=CET, solution=gsModel, simParam=SP)
PYT = setEBV(pop=PYT, solution=gsModel, simParam=SP)
AYT = setEBV(pop=AYT, solution=gsModel, simParam=SP)
UYT = setEBV(pop=UYT, solution=gsModel, simParam=SP)

for(year in (burninYears+1):nCycles){ # beginning of Genomic Selection
  cat("Advancing breeding with broad-adaptation program
      strategy year:",year,"of", nCycles, "\n")

  # Selecting variety for release
  variety <- selectInd(UYT,nInd=nVarietySel, use="ebv",simParam=SP)

  # Advance year of breeding for genomic selection
  # Year 5 Uniform yield trial  (UYT)
  UYT <- selectInd(pop=AYT,  nInd=nUYT, use="ebv", simParam=SP)

  # Invoke the function to phenotype selected UYT clones in 10 locations
  UYT <- gxeSim(pval1=0.1, pval2=0.9, pop=UYT,
                varE=errVarUYT, nreps=repUYT,locSize=8)
  
  UYT <- setEBV(UYT, solution=gsModel, simParam=SP) # estimate GEBV for the CET lines
  
  
  # Advance Yield Trial (AYT)
  AYT = selectInd(pop=PYT, nInd=nAYT,use="ebv",simParam=SP)

  # Invoke the function to phenotype selected AYT clones in 4 locations
  AYT <- gxeSim(pval1=0.3, pval2=0.7,pop=AYT,
                varE=errVarAYT, nreps=repAYT,locSize=4)
  
  AYT <- setEBV(AYT, solution=gsModel) # estimate GEBV for the CET lines

  # Preliminary  Yield  Trial (PYT)
  PYT = selectInd(pop=CET, nInd=nPYT, use="ebv",simParam=SP)

  # Invoke the function to phenotype selected PYT clones in 2 locations
  PYT <- gxeSim(pval1=0.3, pval2=0.7, pop=PYT,
                varE=errVarPYT, nreps=repPYT, locSize=2)
  
  PYT <- setEBV(PYT, solution=gsModel, simParam=SP) # estimate GEBV for the CET lines


  # Clonal Evaluation Trial (CET) - Implement GS
  # Select individuals from the SDN stage based on phenotype performance

  CET <- selectWithinFam(pop=SDN, nInd=famSize, use="pheno",simParam=SP)
  CET <- selectInd(pop=CET, nInd=nCET, simParam=SP)
  CET <- setPheno(pop=CET, varE=errVarCET, reps=repCET,
                  p=0.5, simParam=SP)

  CET <- setEBV(CET, solution=gsModel, simParam=SP) # estimate GEBV for the CET lines

  # Seedling Nursery
  SDN <- setPheno(pop=F1, varE=errVarSDN, reps=repSDN, p=0.5, simParam=SP)

  # recycle new Parents
  # We will use gebv  to select parents among
  # AYT, PYT and CET to potentially boost selection accuracy and reduce
  # generation interval

  # select new parents (50) based on gebv from CET
  parSelCand <- c(UYT,AYT,PYT,CET)
  parSelCand <- setEBV(pop=parSelCand, solution=gsModel,simParam=SP)
  parents <- selectInd(pop=parSelCand, nInd=nInd(parents), use="ebv",simParam=SP)
  F1 <- randCross(pop=parents, nCrosses=nCrosses, nProgeny=nProgeny,
                  simParam=SP)

  # report the result
  meanGV[year] <- meanG(PYT)
  varGV[year] <- varG(PYT)
  
  # Selection accuracy at each evaluation stage
  accCET[year] <-  cor(gv(CET), ebv(CET))
  accPYT[year] <-  cor(gv(PYT), ebv(PYT))
  accAYT[year] <-  cor(gv(AYT), ebv(AYT))
  accUYT[year] <-  cor(gv(UYT), ebv(UYT))
  

  # Update training population and genomic prediction model
  # by retaining the last 2 year of training data
  trainPop <-  c(trainPop[-(1:nInd(c(CET, PYT,AYT,UYT)))],c(CET,PYT,AYT,UYT))
  gsModel <- RRBLUP(pop=trainPop, traits=1, use="pheno",simParam=SP)
}

gGain <- meanGV - meanGV[burninYears]
relVar <- varGV/varGV[burninYears]

simParms <- data.frame(simRun=rep(REP,nCycles),
                     year=(1:nCycles),
                     scenario= rep("Broad adaptation", nCycles),
                     meanGV=meanGV,
                     varGV=varGV,
                     gGain=gGain,
                     relVar=relVar,
                     accCET=accCET,
                     accPYT=accPYT,
                     accAYT=accAYT,
                     accUYT=accUYT)

#write.csv(simParms,file = paste0("./data/broadGxE_parms","_",REP,".csv"),row.names = FALSE)

#saveRDS(result,file = paste0("simulated_data/PYT_GS","_",REP,".rds"))
