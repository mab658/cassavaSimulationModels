# This script implemented conventional breeding scheme bssed on phenotypic selection

for(year in (burninYears+1):nCycles){
  # Desired number of years of a future breeding program
  cat("\n Advance breeding with conventional breeding strategy year:",year,"of", burninYears+futureYears,"\n")
  
  # select the best variety for release
  variety <- selectInd(UYT,nInd=nVarietySel, use="pheno",simParam=SP)

  # Uniform Yield Trial (UYT)
  UYT <- selectInd(pop=AYT, nInd=nUYT, use="pheno",simParam=SP)

  # Invoke the function to phenotype selected UYT clones in 10 locations
  UYT <- gxeSim(pval1=0.1, pval2=0.9, pop=UYT,
                varE=errVarUYT, nreps=repUYT,
                locSize=8)

  # Advance Yield Trial (AYT)
  AYT <-  selectInd(pop=PYT, nInd=nAYT, use="pheno",simParam=SP)

  # Invoke the function to phenotype selected AYT clones in 4 locations
  AYT <- gxeSim(pval1=0.3, pval2=0.7, pop=AYT,
                varE=errVarAYT, nreps=repAYT,locSize=4)


  # Preliminary  Yield  Trial (PYT)
  PYT <- selectInd(pop=CET, nInd=nPYT, use="pheno",simParam=SP)

  # Invoke the function to phenotype selected PYT clones in 2 locations
  PYT <- gxeSim(pval1=0.2, pval2=0.7, pop=PYT,
                varE=errVarPYT,nreps=repPYT,locSize=2)

  # Clonal Evaluation Trial (CET)
  CET <- selectWithinFam(pop=SDN, nInd=famSize, use="pheno", simParam=SP)
  CET <- selectInd(pop=CET, nInd=nCET, use="pheno", simParam=SP)
  CET <- setPheno(pop=CET, varE=errVarCET, reps=repCET,p=0.5, simParam=SP)
  

  # Seedling Nursery (SDN)
  SDN <- setPheno(pop=F1, varE=errVarSDN,reps=repSDN,p=0.5,simParam=SP)


  # Update and recycle parents in the crossing block
  parents <- selectInd(pop=c(UYT,AYT),nInd=50, 
                       use="pheno",simParam=SP)

  # Crossing block - Make crosses to generate new population
  F1 <- randCross(pop=parents, nCrosses=nCrosses,
                  nProgeny=nProgeny, simParam=SP)


  # compute the mean and variance of genetic values over years
  meanGV[year] <- meanG(PYT)
  varGV[year] <- varG(PYT)
  
  # Selection accuracy at each evaluation stage
  accCET[year] <-  cor(gv(CET), pheno(CET))
  accPYT[year] <-  cor(gv(PYT), pheno(PYT))
  accAYT[year] <-  cor(gv(AYT), pheno(AYT))
  accUYT[year] <-  cor(gv(UYT), pheno(UYT))
}

# center mean of GV  and scale variance GV 
# based on final burn-in phase
gGain <- meanGV - meanGV[burninYears]
relVar <- varGV/varGV[burninYears]

simParms <- data.frame(simRun=rep(REP,nCycles),
                       year=(1:nCycles),
                       scenario= rep("Conv", nCycles),
                       meanGV=meanGV,
                       varGV=varGV,
                       gGain=gGain,
                       relVar=relVar,
                       accCET=accCET,
                       accPYT=accPYT,
                       accAYT=accAYT,
                       accUYT=accUYT)

write.csv(simParms,file = paste0("./data/conv","_",REP,".csv"),row.names = FALSE)
#saveRDS(output,file = paste0("simulated_data/Base","_",REP,".rds"))
