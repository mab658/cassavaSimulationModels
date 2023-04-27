# function to fit Variance component estimate model BLUP model
# To estimate broad-sense heritability at UYT stage in conventional breeding scheme

varModel <- function(METdat){
# Model fit to estimate variance component for heritability estimate

  nEnv <- length(unique(METdat$env))
 
  METdat$gen <- as.factor(METdat$gen)
  METdat$env <- as.factor(METdat$env)
 
  modelFit <- asreml(fixed = pheno~env,
                        random =~ gen,
                        residual =~idv(units),
                        na.action=na.method(y='include',x='include'),
                        workspace = 210e07, maxit = 80,data = METdat,trace = F)

  # Loop to ensure model converges
  while (!modelFit$converge) {
        modelFit <- update.asreml(modelFit)
        print("Using update.asreml")
       } # end while loop

  varEst <- summary(modelFit, coef=TRUE)$varcomp # variance component estimates
  varGE <- varEst[2,1] - METdat$errVar[1] # extract out variance of GxE

  H2 <- varEst[1,1]/(varEst[1,1] + varGE/nEnv + METdat$errVar[1]) # compute heritability
  return(H2)
} # end of function
