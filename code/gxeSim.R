# function to simulate GxE at UYT trials
gxeSim <- function(pval1,pval2,genPop,varE,nreps,scenario,locSize){

  GV = matrix(0, nrow=nInd(genPop), ncol=locSize)
  dimnames(GV)<- list(genPop@id,paste(paste0("L",1:locSize)))

  # A p-value of the environmental covariate
  P <- seq(pval1, pval2, length.out=locSize)

  for (loc in 1:locSize){  # begin loop
    GV[,loc] = setPheno(pop=genPop, varE=varE, reps=nreps,
                        p=P[loc],onlyPheno=TRUE, simParam=SP)
  } #end loop

  genPop@pheno <- as.matrix(rowMeans(GV))

  return(genPop)
}
