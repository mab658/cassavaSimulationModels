# function to simulate GxE at UYT trials
gxeSim <- function(pval1,pval2,pop,varE,nreps,locSize){

  GV = matrix(0, nrow=nInd(pop), ncol=locSize)
  dimnames(GV)<- list(pop@id,paste(paste0("L",1:locSize)))

  # A P-value of the environmental covariate
  P <- round(seq(pval1, pval2, length.out=locSize),1)

  for (loc in 1:locSize){  # begin loop
    GV[,loc] = setPheno(pop=pop, varE=varE, reps=nreps,
                        p=P[loc], fixEff=year, onlyPheno=TRUE, simParam=SP)
  } #end loop

  pop@pheno <- as.matrix(rowMeans(GV))

  return(pop)
}
