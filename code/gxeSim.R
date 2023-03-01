# function to simulate phenotypic data for GxE
gxeSim <- function(pval1,pval2,pop,varE,nreps,nLocs){
  stage <- deparse(substitute(pop))

  GV = matrix(0, nrow=nInd(pop), ncol=nLocs)
  dimnames(GV)<- list(pop@id,paste(paste0("L",1:nLocs)))

  # A P-value of the environmental covariate
  P <- round(seq(pval1, pval2, length.out=nLocs),1)

  for (loc in 1:nLocs){  # begin loop
    GV[,loc] = setPheno(pop=pop, varE=varE, reps=nreps, fixEff=year,
                        p=P[loc], onlyPheno=TRUE, simParam=SP)
  } #end loop

  pop@pheno <- as.matrix(rowMeans(GV))
  entryMeanVarE <- varE/nreps/nLocs
  trialRec <- data.frame(id=pop@id,year=year,
                         stage=stage,
                         pheno=pop@pheno,errVar=entryMeanVarE)
  rownames(trialRec) <- NULL
  gxeOut <- list(pop,trialRec)
  return(gxeOut)
}

gxeSimNew <- function(pvalVec,pop,varE,nreps,nLocs){
  # length of pvalVec == nLocs
  stage <- deparse(substitute(pop))

  GV = matrix(0, nrow=nInd(pop), ncol=nLocs)
  dimnames(GV)<- list(pop@id,paste(paste0("L",1:nLocs)))

  for (loc in 1:nLocs){  # begin loop
    GV[,loc] = setPheno(pop=pop, varE=varE, reps=nreps, fixEff=year,
                        p=pvalVec[loc], onlyPheno=TRUE, simParam=SP)
  } #end loop

  pop@pheno <- as.matrix(rowMeans(GV))
  entryMeanVarE <- varE/nreps/nLocs
  trialRec <- data.frame(id=pop@id,year=year,
                         stage=stage,
                         pheno=pop@pheno,errVar=entryMeanVarE)
  rownames(trialRec) <- NULL
  gxeOut <- list(pop,trialRec)
  return(gxeOut)
}
