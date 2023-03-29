# function to simulate phenotypic data for GxE
gxeSim <- function(pvalVec,pop,varE,nreps,locID,year){ # not nLocs
  # length of pvalVec == nLocs
  stage <- deparse(substitute(pop))

  metPheno = matrix(0, nrow=nInd(pop), ncol=length(locID))
  dimnames(metPheno)<- list(pop@id,paste(paste0("L",locID)))


  # A P-value of the environmental covariate

  for (loc in 1:length(locID)){  # begin loop
    metPheno[,loc] = setPheno(pop=pop, varE=varE, reps=nreps,
                        p=pvalVec[loc], onlyPheno=TRUE, simParam=SP) 
  } #end loop

  pop@pheno <- as.matrix(rowMeans(metPheno))
  entryMeanVarE <- varE/nreps/length(locID)
  trialRec <- data.frame(gen=pop@id,year=year,
                         stage=stage,
                         metPheno,errVar=entryMeanVarE)

  rownames(trialRec) <- NULL
 
  trialRec  <- trialRec  %>%
    gather(key = loc,value = pheno,-c(gen, year, stage,errVar))
 
  # concatenate loc and year as environment
  trialRec$env <- paste0(trialRec$loc,"_",trialRec$year)
  trialRec  <- trialRec[,c("loc","year","env","stage","gen","pheno","errVar")] # rearrange the columns
  
  gxePheno <- list(pop,trialRec)
  return(gxePheno)
}
