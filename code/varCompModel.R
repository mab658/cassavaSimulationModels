# function to fit Variance component estimate model BLUP model


varModel <- function(METdat){
	# Model fit to estimate variance component for heritability estimate

  	#nEnv <- length(unique(METdat$env))

	phenoDat  <- METdat  %>%
                group_by(env,gen,stage) %>%
                summarise(pheno=mean(pheno,na.rm=T),errVar=mean(errVar,na.rm=T),.groups = 'drop')

 	phenoDat$gen <- as.factor(phenoDat$gen)
	phenoDat$env <- as.factor(phenoDat$env)
	phenoDat$wt <- 1/phenoDat$errVar
  	modelFit <- asreml(fixed = pheno~env,
                        random =~ gen,
                        residual =~idv(units),
			weights = wt,
                        na.action=na.method(y='include',x='include'),
                        workspace = 210e07, maxit = 80,data = phenoDat,trace = F)

  	# Loop to ensure model converges
  	while (!modelFit$converge) {
        	modelFit <- update.asreml(modelFit)
        	print("Using update.asreml")
       	} # end while loop

  	varEst <- summary(modelFit, coef=TRUE)$varcomp # variance component estimates
  	#varGE <- varEst[2,1] - METdat$errVar[1] # extract out variance of GxE

  	#H2 <- varEst[1,1]/(varEst[1,1] + varGE/nEnv + METdat$errVar[1]) # compute heritability
	H2 <- varEst[1,1]/(varEst[1,1]+varEst[2,1]) # Broad-sense heritability

  	return(H2)
} # end of function
