# function to fit GBLUP (estimate GEBV)  and BLUP (estimate GETGV)  model using inverse of G-matrix in sparse form

gsModel <- function(snpsMarker,pheno){

	#  quality control SNPs marker data
  	# remove duplicates individuals from genotyping data
  	snpsMarker <- snpsMarker[!duplicated(rownames(snpsMarker)),]
  	snpsFilter <- qc.filtering(M=snpsMarker, base=FALSE, ref=NULL,
        marker.callrate=0.2,ind.callrate=0.2, maf=0.05,heterozygosity=0.95,
        Fis=1, impute=FALSE,  plots=FALSE,message=FALSE)$M.clean

  	# compute genomic relationship matrix from filtered SNPS marker
  	Gmat  <- G.matrix(M = snpsFilter, method = "VanRaden", na.string = NA,message=FALSE)$G

 	# Tuning the G-matrix (Gmat) by bending i.e. adjust it near a positive definite
  	# by making some of its very small or negative eigenvalues slightly positive (Matrix::nearPD()).
  	G_bend <- G.tuneup(G=Gmat, bend=TRUE, eig.tol = 1e-06,message=FALSE)$Gb # Bend the Gmat matrix

  	# Tuning the G-matrix (Gmat) by blending # 
  	G_blend <- G.tuneup(G=G_bend, blend=TRUE, pblend=0.02,message=FALSE)$Gb

  	phenoMean  <- pheno  %>%
        	group_by(env,gen,stage) %>%
        	summarise(pheno=mean(pheno,na.rm=T),errVar=mean(errVar,na.rm=T),.groups = 'drop')

  	# match individuals phenotyped to genotyped and vice versa
   	phenoDat <- phenoMean[phenoMean$gen %in% rownames(G_blend),]
   	#G_blend <- G_blend[rownames(G_blend) %in% phenoDat$gen,]
  
	phenoDat$env <- as.factor(phenoDat$env)
  	phenoDat$gen <- as.factor(phenoDat$gen) # coerce genotype to a factor
  	phenoDat$wt <- 1/phenoDat$errVar

  	# get the inverse of blend genomic relationship matrix
  	kinv <<- G.inverse(G = G_blend, sparseform = TRUE,message = FALSE)$Ginv

  	# Train genomic selection (GS) model and safe the model fit
  	modelFit <- asreml(fixed = pheno~env,
     		random =~vm(gen,kinv), # model only additive effect through GRM
        	residual =~idv(units),
                weights = wt,
		na.action=na.method(y='include', x='include'),
                workspace = 210e07, maxit = 80,data = phenoDat,trace = F)

  	# Loop to ensure model converges
  	while (!modelFit$converge) {
        	modelFit <- update.asreml(modelFit)
        	cat("Using update.asreml","\n")
  	}


	# estimate of additive  variance component
	varEst <- summary(modelFit, vparameters=TRUE)$varcomp
	
	h2 <- varEst[1,1]/(varEst[1,1]+varEst[2,1]) # narrow-sense heritability
  	# extract the genomic estimated breeding value (GEBV)
  	genEff <- summary(modelFit, coef=TRUE)$coef.random

	gebv <- genEff[grep("^vm",rownames(genEff)),1] # GEBV i.e additive genetic effect
	names(gebv) <- as.numeric(substr(rownames(genEff),start=15,stop=nchar(rownames(genEff))))
	

	# Fit linear mixed model to estimate total genetic effect or value
        modelFit <- asreml(fixed = pheno~1,
                random =~ idv(gen), # model total genetic effect                        
                residual =~idv(units),
                weights = wt,
                na.action=na.method(y='include', x='include'),
                workspace = 210e07, maxit = 80,data = phenoDat,trace = F)

        # Loop to ensure model converges
        while (!modelFit$converge) {
                modelFit <- update.asreml(modelFit)
                print("Using update.asreml")
        }


	# estimate of additive  variance component
        varEst <- summary(modelFit, vparameters=TRUE)$varcomp

	H2 <- varEst[1,1]/(varEst[1,1]+varEst[2,1]) # Broad-sense heritability

	# extract out the  total genetic value (tgv) from the linear mixed model
        genEff <- summary(modelFit, coef=TRUE)$coef.random
	egv <- genEff[grep("^gen",rownames(genEff)),1] # non additive estimated  genetic effect (egv)
	names(egv) <- as.numeric(substr(rownames(genEff),start=5,stop=nchar(rownames(genEff))))

	return(list(gebv, h2, egv,H2))
} # end of function
