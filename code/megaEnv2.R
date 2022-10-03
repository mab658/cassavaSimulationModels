for(year in (burninYears+1):nCycles){# nCycles=burninYears+futureYears

  cat("Advancing narrow-adaptation Mega environment 2 breeding program
       year:",year,"of", nCycles, "\n")

  # Selecting variety for release
  variety <- selectInd(pop=UYT,nInd=nVarietySel, use="pheno")

  # Year 5 Uniform yield trial  (UYT)
  UYT <- selectInd(pop=AYT,  nInd=nUYT, use="pheno", simParam=SP)

  # Invoke the function to phenotype selected UYT clones in 5 locations
  UYT <- gxeSim(pval1=0.55,pval2=0.95,genPop=UYT,
                varE=errVarUYT, nreps=repUYT,
                scenario="Narrow adaptation",locSize=5)


  # Year 4 Advance Yield Trial (AYT)
  AYT = selectInd(pop=PYT, nInd=nAYT,use="pheno")

  # Invoke the function to phenotype selected AYT clones in 2 locations
  AYT <- gxeSim(pval1=0.55, pval2=0.65, genPop=AYT,
                varE=errVarAYT, nreps=repAYT,
             scenario="Narrow adaptation", locSize=2)

  # Year 3 Preliminary  Yield  Trial (PYT)

  PYT <- selectInd(pop=CET, nInd=nPYT, use="ebv")

  # Invoke the function to phenotype selected PYT clones in 1 locations
  PYT <- gxeSim(pval1=0.55, pval2=0.55, genPop=PYT,
                varE=errVarPYT, nreps=repPYT,
                scenario="Narrow adaptation", locSize=1)

  # Year 2  Clonal evaluation
  # Select individuals from the SDN stage based on phenotype performance
  CET <- selectWithinFam(pop=SDN, nInd=famSize, use="pheno")
  CET <- selectInd(pop=CET, nInd=500, simParam=SP)

  # Invoke the function to phenotype selected PYT clones in 1 locations
  CET <- gxeSim(pval1=0.55, pval2=0.55, genPop=CET,
                varE=errVarCET, nreps=repCET,
             scenario="Narrow adaptation", locSize=1)


  CET = setEBV(CET, gsModel) # estimate breeding values for the CET lines

  # Evaluate accuracy of selection based on phenotype values

  accPheno[year] = cor(gv(CET), pheno(CET))
  # Evaluate accuracy of selection based on genomic EBV
  accEbv[year] = cor(gv(CET), ebv(CET))

  # Year 1  seedling nursery
  SDN <- F1

  SDN <- setPheno(pop=SDN, varE=errVarSDN, reps=repSDN, p=0.55,
                  fixEff=year, simParam=SP)

  # select new parents for crossing (50)
  parents <- selectInd(pop=c(UYT,AYT), nInd=nInd(parents), use="ebv")
  F1 <- randCross(pop=parents, nCrosses=100, nProgeny=nProgeny,
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
                       scenario= rep("Narrow adaptation", nCycles),
                       meanGV=meanGV,
                       varGV=varGV,
                       gGain=gGain,
                       relVar=relVar,
                       accuracy=accEbv)

write.csv(simParms,file=paste0("../data/narrowME2_parms","_",REP,".csv"),row.names=FALSE)
#saveRDS(result,file = paste0("simulated_data/PYT_GS","_",REP,".rds"))

