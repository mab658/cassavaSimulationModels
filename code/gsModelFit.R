# function to fit GBLUP (estimate GEBV)  and BLUP (estimate GETGV)  model using inverse of G-matrix in sparse form

gsModel <- function(snpsMarker,pheno){

	#  quality control SNPs marker data
  	# remove duplicates individuals from genotyping data
  	snpsMarker <- snpsMarker[!duplicated(rownames(snpsMarker)),]
  	
	snpsFilter  <- qc.filtering(M=snpsMarker, base=FALSE, ref=NULL,
        marker.callrate=0.2,ind.callrate=0.2, maf=0.05,heterozygosity=0.95,
        Fis=1, impute=FALSE,  plots=FALSE,message=FALSE)$M.clean

  	# compute genomic relationship matrix from filtered SNPS marker
  	Gmat  <- G.matrix(M = snpsFilter, method = "VanRaden", na.string = NA,message=FALSE)$G

 	# Tuning the G-matrix (Gmat) by bending i.e. adjust it near a positive definite
  	# by making some of its very small or negative eigenvalues slightly positive (Matrix::nearPD()).
  	G_bend <- G.tuneup(G=Gmat, bend=TRUE, eig.tol = 1e-06,message=FALSE)$Gb # Bend the Gmat matrix

  	# Tuning the G-matrix (Gmat) by blending # 
  	G_blend <- G.tuneup(G=G_bend, blend=TRUE, pblend=0.02,message=FALSE)$Gb

	# get the inverse of blend genomic relationship matrix
        kinv <<- G.inverse(G = G_blend, sparseform = TRUE,message = FALSE)$Ginv
	

	# order the phenotypic data for modeling
	phenoDat  <- pheno  %>%
		select(env, gen, stage,pheno, errVar) %>%
                arrange(env,gen,stage) 

	phenoDat$env <- as.factor(phenoDat$env)
  	phenoDat$gen <- as.factor(phenoDat$gen) # coerce genotype to a factor
  	phenoDat$wt <- 1/phenoDat$errVar

  	# Train genomic selection (GS) model and safe the model fit
  	modelFit <- asreml(fixed = pheno~env,
     		random =~vm(gen,kinv) + ide(gen), # model only additive effect through GRM
        	residual =~idv(units),
                weights = wt,
		na.action=na.method(x='include'),
                workspace = 210e07, maxit = 100,data = phenoDat,trace = F)

  	# Loop to ensure model converges
  	while (!modelFit$converge) {
        	modelFit <- update.asreml(modelFit)
        	cat("Using update.asreml","\n")
  	}


	# estimate of additive  variance component
	varEst <- summary(modelFit, vparameters=TRUE)$varcomp
	
	h2 <- varEst[1,1]/(varEst[1,1]+varEst[2,1]+varEst[3,1]) # narrow-sense heritability
	H2 <- (varEst[1,1] + varEst[2,1])/(varEst[1,1]+varEst[2,1]+varEst[3,1]) # Broad-sense heritability

  	# extract the genomic estimated breeding value (GEBV)
  	genEff <- summary(modelFit, coef=TRUE)$coef.random

	gebv <- genEff[grep("^vm",rownames(genEff)),1] # GEBV i.e additive genetic effect
	names(gebv) <- as.numeric(substr(names(gebv),start=15,stop=nchar(names(gebv))))

	# non-additive effect
	nonAddEff <- genEff[grep("^ide",rownames(genEff)),1] #  i.e nonadditive genetic effect
        names(nonAddEff) <- as.numeric(substr(names(nonAddEff),start=10,stop=nchar(names(nonAddEff))))
	
	egv <- gebv + nonAddEff # sum additive and non additive effects
	return(list(gebv, h2, egv,H2))
} # end of function
