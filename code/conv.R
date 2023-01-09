# This script implemented conventional breeding scheme bssed on phenotypic selection

for(year in (burninYears+1):nCycles){
  # Desired number of years of a future breeding program
  cat("\n Advance breeding with conventional breeding strategy year:",year,"of", burninYears+futureYears,"\n")
  
  # select the best variety for release

  # Selection accuracy
  accUYT[year] <-  cor(gv(UYT), pheno(UYT))
  variety <- selectInd(pop=UYT,nInd=nVarietySel, use="pheno",simParam=SP)

  
  # Uniform Yield Trial (UYT)

  # Selection accuracy
  accAYT[year] <-  cor(gv(AYT), pheno(AYT))
  UYT <- selectInd(pop=AYT, nInd=nUYT, use="pheno",simParam=SP)
  
  # Invoke the function to phenotype selected UYT clones in 10 locations
  UYTrec <- gxeSim(pval1=0.1, pval2=0.9, pop=UYT,
                varE=errVarUYT, nreps=repUYT,
                nLocs=8)

  UYT <- UYTrec[[1]]

  # Selection accuracy
  accPYT[year] <-  cor(gv(PYT), pheno(PYT))

  # Advance Yield Trial (AYT)
  AYT <-  selectInd(pop=PYT, nInd=nAYT, use="pheno",simParam=SP)

  # Invoke the function to phenotype selected AYT clones in 4 locations
  AYTrec <- gxeSim(pval1=0.3, pval2=0.7, pop=AYT,
                varE=errVarAYT, nreps=repAYT,nLocs=4)
  
  AYT <- AYTrec[[1]]

  # Selection accuracy
   accCET[year] <-  cor(gv(CET), pheno(CET))
  
  # Preliminary  Yield  Trial (PYT)
  PYT <- selectInd(pop=CET, nInd=nPYT, use="pheno",simParam=SP)

  # Invoke the function to phenotype selected PYT clones in 2 locations
  PYTrec <- gxeSim(pval1=0.3, pval2=0.7, pop=PYT,
                varE=errVarPYT,nreps=repPYT,nLocs=2)


  PYT <- PYTrec[[1]]

  # Clonal Evaluation Trial (CET)
  CET <- selectWithinFam(pop=SDN, nInd=famSize, use="pheno", simParam=SP)
  CETrec <- gxeSim(pval1=0.5, pval2=0.5, pop=CET,
                varE=errVarCET,nreps=repCET,nLocs=1)

  CET <- CETrec[[1]]


  # Seedling Nursery (SDN)
  SDN <- setPheno(pop=F1, varE=errVarSDN,reps=repSDN,p=0.5,simParam=SP)

  # Update and recycle parents in the crossing block
  parents <- selectInd(pop=c(UYT,AYT),nInd=50,
                       use="pheno",simParam=SP)

  # Crossing block - Make crosses to generate new population
  F1 <- randCross(pop=parents, nCrosses=nCrosses,
                  nProgeny=nProgeny, simParam=SP)


  # save mean & variance of genetic values to track progress of breeding programs
  meanGV[year] <- meanG(PYT)
  varGV[year] <- varG(PYT)  
} # end loop

# put the simulated parameters as a data frame object
simParms <- data.frame(simRun=rep(REP,nCycles),
                       year=(1:nCycles),
                       scenario= rep("Conv", nCycles),
                       meanGV=meanGV,
                       varGV=varGV,
		       varGxE=varGE,

                       accCET=accCET,
                       accPYT=accPYT,
                       accAYT=accAYT,
                       accUYT=accUYT)
