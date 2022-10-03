# Alternative breeding scenario one
# This script implemented broad adaptation  breeding programs
# with genomic selection at CET to boost the selection accuracy
# without changing the generation interval


# Train genomic selection (GS) model and safe the model fit
gsModel <- RRBLUP(pop=trainPop, traits=1, simParam=SP)

# Estimate breeding values for CET, PYT, and AYT for being used to
# constitute the training population set and they are genotyped
CET = setEBV(pop=CET, solution=gsModel)
PYT = setEBV(pop=PYT, solution=gsModel)
AYT = setEBV(pop=AYT, solution=gsModel)
UYT = setEBV(pop=UYT, solution=gsModel)

for(year in (burninYears+1):nCycles){# nCycles=burninYears+futureYears

  cat("Advancing breeding with broad-adaptation program
      strategy year:",year,"of", nCycles, "\n")

  # Selecting variety for release
  variety <- selectInd(UYT,nInd=nVarietySel, use="pheno")

  # Advance year of breeding for genomic selection
  # Year 5 Uniform yield trial  (UYT)
  UYT <- selectInd(pop=AYT,  nInd=nUYT, use="pheno", simParam=SP)

  # Invoke the function to phenotype selected UYT clones in 10 locations
  UYT <- gxeSim(pval1=0.05, pval2=0.95, genPop=UYT,
                varE=errVarUYT, nreps=repUYT,
                scenario="Broad adaptation",locSize=10)


  # Year 4 Advance Yield Trial (AYT)
  AYT = selectInd(pop=PYT, nInd=nAYT,use="pheno")

  # Invoke the function to phenotype selected AYT clones in 4 locations
  AYT <- gxeSim(pval1=0.35, pval2=0.65, genPop=AYT,
                varE=errVarAYT, nreps=repAYT,
                scenario="Broad adaptation", locSize=4)


  # select the best PYT lines based on estimated breeding values (EBV)
  # from the phenotyped and genotyped evaluated CET lines

  # Year 3 Preliminary  Yield  Trial (PYT)
  PYT = selectInd(pop=CET, nInd=nPYT, use="ebv")

  # Invoke the function to phenotype selected PYT clones in 2 locations
  PYT <- gxeSim(pval1=0.45, pval2=0.55, genPop=PYT,
                varE=errVarPYT, nreps=repPYT,
             scenario="Broad adaptation", locSize=2)


  # Year 2 Clonal Evaluation Trial (CET) - Implement GS

  # Select individuals from the SDN stage based on phenotype performance

  CET <- selectWithinFam(pop=SDN, nInd=famSize, use="pheno")
  CET <- selectInd(pop=CET, nInd=nCET, simParam=SP)

  # Invoke the function to phenotype selected CET clones in 1 locations
  CET <- gxeSim(pval1=0.45, pval2=0.45, genPop=CET,
                varE=errVarCET, nreps=repCET,
             scenario="Broad adaptation", locSize=1)

  CET = setEBV(CET, gsModel) # estimate breeding values for the CET lines

  # Evaluate accuracy of selection based on phenotype values
  accPheno[year] = cor(gv(PYT), pheno(PYT))
  # Evaluate accuracy of selection based on genomic EBV
  accEbv[year] = cor(gv(PYT), ebv(PYT))


  # Year 1  Seedling Nursery
  SDN <- F1
  SDN <- setPheno(pop=SDN, varE=errVarSDN, reps=repSDN, p=0.45,
                  fixEff=year, simParam=SP)

  # recycle new Parents
  # We will use gebv  to select parents among
  # AYT, PYT and CET to potentially boost selection accuracy and reduce
  # generation interval

  # select new parents (50) based on gebv from CET
  parents <- selectInd(pop=c(UYT,AYT), nInd=nInd(parents), use="ebv")
  F1 <- randCross(pop=parents, nCrosses=nCrosses, nProgeny=nProgeny,
                  simParam=SP)

  # report the result
  meanGV[year] <- meanG(PYT)
  varGV[year] <- varG(PYT)

  # Update training population and genomic prediction model
  # by retaining the last 2 year of training data
  trainPop <-  c(trainPop[-(1:nInd(c(CET, PYT,AYT,UYT)))],c(CET,PYT,AYT,UYT))
  gsModel <- RRBLUP(pop=trainPop, traits=1, simParam=SP)
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
                     accuracy=accEbv)

write.csv(simParms,file = paste0("../data/broadGxE_parms","_",REP,".csv"),row.names = FALSE)

#saveRDS(result,file = paste0("simulated_data/PYT_GS","_",REP,".rds"))
