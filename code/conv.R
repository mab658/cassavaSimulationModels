# This script implemented conventional breeding scheme
# without genomic selection

for(year in (burninYears+1):nCycles){
  # Desired number of years of a future breeding program
  cat("\n Advance breeding with conventional breeding strategy year:",year,"of", burninYears+futureYears,"\n")

  variety <- selectInd(UYT,nInd=nVarietySel, use="pheno")

  # Year 5 Uniform Yield Trial (UYT)
  UYT = selectInd(pop=AYT, nInd=nUYT, use="pheno")

  # Invoke the function to phenotype selected UYT clones in 10 locations
  UYT <- gxeSim(pval1=0.05, pval2=0.95, genPop=UYT,
                varE=errVarUYT, nreps=repUYT,
                scenario="Conv", locSize=10)

  # Year 4 Advance Yield Trial (AYT)
  AYT = selectInd(pop=PYT, nInd=nAYT, use="pheno")

  # Invoke the function to phenotype selected AYT clones in 4 locations
  AYT <- gxeSim(pval1=0.35, pval2=0.65, genPop=AYT,
                varE=errVarAYT, nreps=repAYT,
                scenario="Conv", locSize=4)


  # Year 3 Preliminary  Yield  Trial (PYT)
  PYT = selectInd(pop=CET, nInd = nPYT, use = "pheno")

  # Invoke the function to phenotype selected PYT clones in 2 locations
  PYT <- gxeSim(pval1=0.45, pval2=0.55, genPop=PYT,
                varE=errVarPYT,nreps=repPYT,
             scenario="Conv", locSize=2)

  # Year 2 Clonal Evaluation Trial (CET)
  CET <- selectWithinFam(pop=SDN, nInd=famSize, use="pheno")
  CET <- selectInd(pop=CET, nInd=nCET, use="pheno", simParam=SP)

  # Invoke the function to phenotype selected CET clones in 1 locations
  CET <- gxeSim(pval1=0.45, pval2=0.45, genPop=CET,
                varE=errVarCET, nreps=repPYT,
                scenario="Conv", locSize=1)


  # Year 1  Seedling Nursery
  SDN <- F1
  SDN <- setPheno(pop=SDN, varE=errVarSDN,reps=repSDN, p=0.45,simParam=SP)


  # Update and recycle parents in the crossing block
  parents <- selectInd(pop=c(UYT,AYT),nInd=50, use="pheno")

  # Crossing block - Make crosses to generate new population
  F1 <- randCross(pop=parents, nCrosses=nCrosses,
                  nProgeny=nProgeny, simParam=SP)


  # compute the mean and variance of genetic values over years
  meanGV[year] <- meanG(PYT)
  varGV[year] <- varG(PYT)
  accPheno[year] <- cor(PYT@gv, PYT@pheno) #report accuracies of phenotypic sel.
}

gGain <- meanGV - meanGV[burninYears]
relVar <- varGV/varGV[burninYears]

simParms <- data.frame(simRun=rep(REP,nCycles),
                       year=(1:nCycles),
                       scenario= rep("Conv", nCycles),
                       meanGV=meanGV,
                       varGV=varGV,
                       gGain=gGain,
                       relVar=relVar,
                       accuracy=accPheno)

write.csv(simParms,file = paste0("../data/conv","_",REP,".csv"),row.names = FALSE)
#saveRDS(output,file = paste0("simulated_data/Base","_",REP,".rds"))
