# function to fit BLUP model, estimate variance component, heritability and predictive ability

varModel <- function(METdat){

	# order the phenotypic data for modeling
        phenoDat  <- METdat  %>%
                select(env, gen, stage,pheno, errVar) %>%
                arrange(env,gen,stage) 

	phenoDat$env <- as.factor(phenoDat$env)
 	phenoDat$gen <- as.factor(phenoDat$gen)
	phenoDat$wt <- 1/phenoDat$errVar
  	modelFit <- asreml(fixed = pheno~env,
                        random =~ gen,
                        residual =~idv(units),
			weights = wt,
                        na.action=na.method(x='include'),
                        workspace = 210e07, maxit = 100,data = phenoDat,trace = F)

  	# Loop to ensure model converges
  	while (!modelFit$converge) {
        	modelFit <- update.asreml(modelFit)
        	print("Using update.asreml")
       	} # end while loop


	# extract variance component estimates and calculate H2
  	varEst <- summary(modelFit, coef=TRUE)$varcomp # variance component estimates
	H2 <- varEst[1,1]/(varEst[1,1]+varEst[2,1]) # Broad-sense heritability

  	return(H2)
} # end of function
