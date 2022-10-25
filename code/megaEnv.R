megaEnv <- function(pval1UYT,pval2UYT,
                    pval1AYT, pval2AYT,
                    pvalPYT,pvalCET,pvalSDN){
  
  for(year in (burninYears+1):nCycles){
    variety <- selectInd(pop=UYT,nInd=nVarietySel, 
                         use="ebv",simParam=SP)
    
    # Uniform Yield Trial (UYT)
    UYT <<- selectInd(pop=AYT,  nInd=nUYT, 
                      use="ebv", simParam=SP)
    
    # Invoke the function to phenotype selected UYT clones in 5 locations
    #pval1=0.05
    UYT <- gxeSim(pval1=pval1UYT,pval2=pval2UYT,pop=UYT,
                  varE=errVarUYT, nreps=repUYT,
                  locSize=5)
    
    UYT <- setEBV(UYT, gsModel)
    
    
    #  Advance Yield Trial (AYT)
    AYT <- selectInd(pop=PYT, nInd=nAYT,use="ebv",
                    simParam=SP)
    
    # Invoke the function to phenotype selected AYT clones in 2 locations
    AYT <- gxeSim(pval1=pval1AYT, pval2=pval2AYT,pop=AYT,
                  varE=errVarAYT, nreps=repAYT,
                  locSize=2)
    
    # estimate breeding values for AYT
    AYT <- setEBV(AYT, gsModel)
    
    
    # Year 3 Preliminary  Yield  Trial (PYT)
    PYT <- selectInd(pop=CET, nInd=nPYT, use="ebv", simParam=SP)
    
    PYT <- setPheno(pop=PYT, varE=errVarPYT, reps=repPYT, 
                    p=pvalPYT, simParam=SP)
    
    # estimate breeding values for CET
    PYT <- setEBV(PYT, gsModel)
    
    # Clonal evaluation trial (CET)
    # Select individuals from the SDN stage based on phenotype performance
    
    CET <- selectWithinFam(pop=SDN, nInd=famSize, use="pheno",simParam=SP)
    CET <- selectInd(pop=CET, nInd=500, use="pheno", simParam=SP)
    
    CET <- setPheno(pop=CET, varE=errVarCET, reps=repCET, 
                    p=pvalCET, simParam=SP)
    
    # estimate breeding values for CET
    CET <- setEBV(CET, gsModel)
    
    
    # seedling nursery
    SDN <- setPheno(pop=F1, varE=errVarSDN, reps=repSDN, 
                    p=pvalSDN,simParam=SP)
    
    # select new parents for crossing (50)
    parents <- selectInd(pop=CET, nInd=nInd(parents), use="ebv",simParam=SP)
    
    F1 <- randCross(pop=parents, nCrosses=100, 
                    nProgeny=nProgeny,
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
    trainPop <-  c(trainPop[-(1:nInd(c(CET,PYT,AYT,UYT)))],c(CET,PYT,AYT,UYT))
    gsModel <- RRBLUP(pop=trainPop, traits=1,use="pheno", simParam=SP)
  } # end loop
  
  gGain <- meanGV - meanGV[burninYears]
  relVar <- varGV/varGV[burninYears]
  
  simParms <- data.frame(simRun=rep(REP,nCycles),
                         year=(1:nCycles),
                         scenario= rep("Narrow adaptation", nCycles),
                         meanGV=meanGV,
                         varGV=varGV,
                         gGain=gGain,
                         relVar=relVar,
                         accCET=accCET,
                         accPYT=accPYT,
                         accAYT=accAYT,
                         accUYT=accUYT)
  return(simParms)
} # end function
