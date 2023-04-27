gsBroad <- function(REP){
 	# Scenario 2
  	# This script implemented broad adaptation breeding programs
  	# with GS to estimate BV in CET based on genomic data to boost 
  	# the selection accuracy at this early stage where phenotype 
  	# information is limited
  	# GS is used to advance from CET, PYT, AYT, and UYT

	
  	cat("\n Modeling scenario 2 broad adapatation enabled GS","\n")

  	cat("\n Simulation run:",REP,"of", nSimRun,"replication", "\n")

  	# load the state of simulation at last year of burnin phase (year=10) for the second scenario
  	load(paste0("./data/burnin_",REP,".rda"))


  	# train the GBLUP model for genomic prediction
   	modelParms <- gsModel(snpsMarker = trainRec$genoBA,pheno = trainRec$phenoBA)


	for(year in (burninYears+1):nCycles){ # beginning of Genomic Selection
  		cat("Replication", REP, "Advancing breeding with broad-adaptation program
      		strategy year:",year,"of", nCycles, "\n")

	        # Use GP model to assign gebv to individuals in the breeding population

		CET@ebv <- as.matrix(modelParms[[1]][CET@id])
        	PYT@ebv <- as.matrix(modelParms[[1]][PYT@id])
        	AYT@ebv <- as.matrix(modelParms[[1]][AYT@id])
        	UYT@ebv <- as.matrix(modelParms[[1]][UYT@id])
		
		# constitute parental selection candidate
		parSelCand <- c(UYT,AYT,PYT,CET)

                if (any(is.na(parSelCand@ebv))){
                        fileName <- paste0("Missing_parentEBV_rep", REP, "_year", year, ".rds")
                        saveRDS(mget(ls()), fileName)
			parSelCand <- parSelCand[!is.na(parSelCand@ebv)] # remove missing value
                }

		parents <- selectInd(pop = parSelCand, nInd = nParents,trait = 1, use = "ebv",simParam = SP)
	
		# number of individuals constitute parental candidates per trial stage
                nParCET[year] <- sum(parents@id %in% CET@id)
                nParPYT[year] <- sum(parents@id %in% PYT@id)
                nParAYT[year] <- sum(parents@id %in% AYT@id)
                nParUYT[year] <- sum(parents@id %in% UYT@id)

		if (any(is.na(CET@ebv))){
      			fileName <- paste0("MissingEBV_CET_rep", REP, "_year", year, ".rds")
      			saveRDS(mget(ls()), fileName)
			CET <- CET[!is.na(CET@ebv)] # remove missing value
  		}


		# Assign GETGV to the population from GP model as surrogates of selection values	
        	PYT@ebv <- as.matrix(modelParms[[3]][PYT@id])
        	AYT@ebv <- as.matrix(modelParms[[3]][AYT@id])
        	UYT@ebv <- as.matrix(modelParms[[3]][UYT@id])


		# assign p-value to each location
  		for (loc in 1:9){
       			allLocPvals[loc] <- runif(1, pvalRange[loc, 1], pvalRange[loc, 2])
  		}  

  		# Selecting from UYT lines for variety for release and its accuracy
		accUYT[year] <-  cor(gv(UYT), ebv(UYT))
  		variety <- selectInd(pop = UYT, nInd = nVarietySel,trait = 1, use = "ebv",simParam = SP)


		# Selecting from AYT lines based on GETGV and its accuracy
		accAYT[year] <-  cor(gv(AYT), ebv(AYT))
  		UYT <- selectInd(pop = AYT, nInd = nUYT, trait = 1, use = "ebv", simParam = SP)


  		# Invoke the function to phenotype selected UYT clones in 8 locations
  
  		UYTrec <- gxeSim(pvalVec=allLocPvals[-5], pop = UYT,
                	varE = errVarUYT, nreps = repUYT, locID=c(1:4,6:9),year=year)
  		UYT <- UYTrec[[1]]  


		# Advance Yield Trial (AYT) - use GETGV to select the best from PYT lines  and its accuracy

		accPYT[year] <-  cor(gv(PYT), ebv(PYT))
  		AYT <- selectInd(pop = PYT, nInd = nAYT,trait = 1,use = "ebv",simParam = SP)
 
  		# Invoke the function to phenotype selected AYT clones in 4 locations
  		AYTrec <- gxeSim(pvalVec=allLocPvals[c(3,4,6,7)],pop = AYT,
                	varE = errVarAYT, nreps = repAYT,  locID = c(3,4,6,7),year=year)
  		AYT <- AYTrec[[1]]


  		# Preliminary  Yield  Trial (PYT) - use GEBV to select the best CET lines and its accurcy
	
		accCET[year] <-  cor(gv(CET), ebv(CET))
  		PYT <- selectInd(pop = CET, nInd = nPYT, use = "ebv", simParam = SP)
 		# Invoke the function to phenotype selected PYT clones in 2 locations
  		PYTrec <- gxeSim(pvalVec=allLocPvals[c(3,7)],pop = PYT,
                	varE = errVarPYT, nreps = repPYT, locID=c(3,7),year=year)
  		PYT <- PYTrec[[1]]  

  
  		# Clonal Evaluation Trial (CET) - advanced individuals from SDN to CET

   		CET <- selectWithinFam(pop = SDN, nInd = famSize, trait = 1, use = "pheno", simParam = SP)
   		CET <- selectInd(pop = CET, nInd = nCET, trait = 1, use = "pheno", simParam = SP)

  		#Invoke the function to phenotype selected SDNclones in 1 locations
  		CETrec <- gxeSim(pvalVec=allLocPvals[5], pop = CET,
                   	varE = errVarCET,nreps = repCET, locID = 5,year=year)  
  		CET <- CETrec[[1]]


  		# Seedling Nursery
  		SDNrec  <- gxeSim(pvalVec=allLocPvals[5], pop = F1,
                  	 varE = errVarSDN,nreps = repSDN, locID = 5,year=year)

  		SDN <- SDNrec[[1]]

 
  		# Update training population
  		# by retaining the last 2 year of training data

  		# Updates the phenotypic data
  		temp <- rbind(CETrec[[2]],PYTrec[[2]], AYTrec[[2]], UYTrec[[2]])

  		# The number 1, 2, 4, and 8 denotes the number of locations per CET, PYT, AYT, and UYT respectively

  		trainRec$phenoBA <- rbind(trainRec$phenoBA[-(1:(nInd(CET)*1 + nInd(PYT)*2 +
                                        nInd(AYT)*4 + nInd(UYT)*8)),],temp)

  		# update genotyping data
  		trainRec$genoBA <- rbind(trainRec$genoBA[-(1:nInd(c(CET, PYT,AYT,UYT))),],
  		pullSnpGeno(pop = c(CET, PYT, AYT, UYT), simParam = SP))

  		# train the GS model for genomic prediction
  		modelParms <- gsModel(snpsMarker = trainRec$genoBA, pheno = trainRec$phenoBA)


  		# crossing block
  		F1 <- randCross(pop = parents, nCrosses = nCrosses, nProgeny = nProgeny, simParam = SP)

  		# to track progress of breeding program at PYT stage
		# Save mean and variance of genetic values
  		meanGV[year] <- meanG(PYT)
 		varGV[year] <- varG(PYT)

		# heritability estimates
       		h2[year] <- modelParms[[2]] # narrow-sense
		H2[year] <- modelParms[[4]] # broad-sense
	} # end loop

	# put the simulated parameters as a data frame object
  	simParms <- data.frame(simRun=rep(REP,nCycles),year=(1:nCycles),scenario= rep("Broad adaptation", nCycles),
		meanGV=meanGV,
     		varGV=varGV,
        	varGxE=varGxE,
		h2=h2,
 		H2 = H2,

 		accCET=accCET,
 		accPYT=accPYT,
     		accAYT=accAYT,
     		accUYT=accUYT,

     		nParCET=nParCET,
     		nParPYT=nParPYT,
     		nParAYT=nParAYT,
		nParUYT=nParUYT)

	return(simParms)

} # end of function
