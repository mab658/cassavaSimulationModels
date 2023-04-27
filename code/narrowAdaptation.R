# A function defined to simulate narrow adaptation
# Test the new breeding scenario for 10-cycle of selection after burnin phase
gsNarrow  <- function(REP,ME,pvalUYT,locUYT,
	            pvalAYT, locAYT,
		    pvalPYT,locPYT){

	cat("\n Narrow breeding program for mega environment ",ME, "\n")

  	cat("\n Simulation run:",REP,"of", nSimRun,"replication", "\n")

  	# load the burnin phase for each replicate of this scenario
  	load(paste0("./data/burnin_",REP,".rda"))

  	if(ME == "ME1"){
        	trainRec$phenoNA <- trainRec$phenoME1
        	trainRec$genoNA <- trainRec$genoME1
		pval=0.36
 	} else if(ME == "ME2"){
        	trainRec$phenoNA <- trainRec$phenoME2
        	trainRec$genoNA <- trainRec$genoME2
   		pval=0.64
  	}


	# invoke gsModel function to fit GBLUP model
	modelParms <- gsModel(snpsMarker = trainRec$genoNA, pheno  = trainRec$phenoNA)


	for(year in (burninYears+1):nCycles){ # begin loop of 10-cycle selection
   

		# Use GP model to assign gebv to individuals in the breeding population
        	CET@ebv <- as.matrix(modelParms[[1]][CET@id])
        	PYT@ebv <- as.matrix(modelParms[[1]][PYT@id])
        	AYT@ebv <- as.matrix(modelParms[[1]][AYT@id])
        	UYT@ebv <- as.matrix(modelParms[[1]][UYT@id])

		# Recycle  new parents in the current year based on gebv before advancing materials
        	parSelCand <- c(UYT,AYT,PYT,CET)

        	if (any(is.na(parSelCand@ebv))){
                	fileName <- paste0("Missing_parentEBV_rep", REP, "_year", year, ".rds")
                	saveRDS(mget(ls()), fileName)
			parSelCand <- parSelCand[!is.na(parSelCand@ebv)] # remove missing value
        	}
        	parents <- selectInd(pop = parSelCand, nInd = nParents,trait = 1, use = "ebv",simParam = SP)
	
		# number of individuals by  trial stages in the parental candidate
        	nParCET[year] <- sum(parents@id %in% CET@id)
        	nParPYT[year] <- sum(parents@id %in% PYT@id)
        	nParAYT[year] <- sum(parents@id %in% AYT@id)
        	nParUYT[year] <- sum(parents@id %in% UYT@id)

		if (any(is.na(CET@ebv))){
                	fileName <- paste0("MissingEBV_CET_rep", REP, "_year", year, ".rds")
                	saveRDS(mget(ls()), fileName)
                	CET <- CET[!is.na(CET@ebv)] # remove missing value
        	}


		# assign p-value to each location
  		for (loc in 1:9){
       			allLocPvals[loc] <- runif(1, pvalRange[loc, 1], pvalRange[loc, 2])
  		}
 
    	
		# Assign GETGV to the population from GP model as surrogates of selection values
                PYT@ebv <- as.matrix(modelParms[[3]][PYT@id])
                AYT@ebv <- as.matrix(modelParms[[3]][AYT@id])
                UYT@ebv <- as.matrix(modelParms[[3]][UYT@id])


		cat("Replication", REP, "Advancing breeding with narrow-adaptation program",ME,
        	"year:",year,"of", nCycles, "\n")

    	
		# Selecting from UYT lines based on GETGV for variety for release and its accuracy

		accUYT[year] <- cor(gvADG(UYT,pval), ebv(UYT))
		#accUYT[year] <-  cor(gv(UYT), ebv(UYT))
    		variety <- selectInd(pop = UYT,nInd = nVarietySel,trait = 1,use = "ebv",simParam = SP)

  
    		# Uniform Yield Trial (UYT) - use GETGV to select the best from AYT lines to UYT
		accAYT[year] <- cor(gvADG(AYT,pval),ebv(AYT))
		#accAYT[year] <-  cor(gv(AYT), ebv(AYT))
    		UYT <- selectInd(pop = AYT, nInd = nUYT, trait = 1, use = "ebv", simParam = SP)

    		# Invoke the function to phenotype selected UYT clones in 4 locations
    		UYTrec <- gxeSim(pvalVec = pvalUYT, pop=UYT,
			varE = errVarUYT, nreps = repUYT,locID = locUYT,year=year)
   
    		UYT <- UYTrec[[1]] # extract out UYT pop from the list output
 
  
  		#  Advance Yield Trial (AYT) - use GETGV to select the best PYT lines and its accuracy
		accPYT[year] <- cor(gvADG(PYT,pval), ebv(PYT))
		#accPYT[year] <-  cor(gv(PYT), ebv(PYT))
    		
		AYT <- selectInd(pop = PYT, nInd = nAYT,trait = 1, use = "ebv", simParam = SP)
    		# Invoke the function to phenotype selected AYT clones in 2 locations
    		AYTrec <- gxeSim(pvalVec = pvalAYT,pop = AYT,
                 	 varE = errVarAYT, nreps = repAYT, locID = locAYT,year=year)
    
    		AYT <- AYTrec[[1]]

    		# Preliminary Yield Trial (PYT) - use GEBV to select the best CET lines to  PYT
		accCET[year] <- cor(gvADG(CET,pval), ebv(CET))    
		#accCET[year] <-  cor(gv(CET), ebv(CET))
    		PYT <- selectInd(pop = CET, nInd = nPYT, trait = 1, use = "ebv", simParam=SP)   
    		PYTrec <- gxeSim(pvalVec = pvalPYT, pop = PYT,
                     varE = errVarPYT,nreps = repPYT,locID = locPYT,year=year)
  

    		# Clonal Evaluation Trial (CET) - Implement GS
    		# Select individuals from the SDN stage based on phenotype performance

    		CET <- selectWithinFam(pop = SDN, nInd = famSize, trait = 1, use = "pheno", simParam = SP)
    		CET <- selectInd(pop = CET, nInd = 500, trait = 1, use = "pheno", simParam = SP)
    		CETrec <- gxeSim(pvalVec = allLocPvals[5], pop = CET,
                     varE=errVarCET,nreps=repCET,locID = 5,year=year)
    
    		CET <- CETrec[[1]] # CET population extract from the list 


    		# seedling nursery
    		SDNrec  <- gxeSim(pvalVec=allLocPvals[5], pop = F1,
                   varE = errVarSDN,nreps = repSDN, locID = 5,year=year)

    		SDN <- SDNrec[[1]]

    		# Update training population
    		# by retaining the last 2 year of training data

		if(year <= (burninYears + limityear)){

        		trainRec$phenoNA <- trainRec$phenoNA[-(1:(nCET + nInd(PYT) + nInd(AYT)*2 + nInd(UYT)*4)),]
        		trainRec$genoNA <- trainRec$genoNA[-(1:(nCET + nInd(c(PYT,AYT,UYT)))),]
        	} else {
        		trainRec$phenoNA <- trainRec$phenoNA[-(1:(nInd(CET) + nInd(PYT) + nInd(AYT)*2 + nInd(UYT)*4)),]
        		trainRec$genoNA <- trainRec$genoNA[-(1:nInd(c(CET, PYT,AYT,UYT))),]
        	}

   		# Updates the phenotyping and genotyping data
    		temp <- rbind(CETrec[[2]],PYTrec[[2]],AYTrec[[2]],UYTrec[[2]])

    		trainRec$phenoNA <- rbind(trainRec$phenoNA,temp)
    		trainRec$genoNA <- rbind(trainRec$genoNA,pullSnpGeno(pop=c(CET, PYT, AYT, UYT),simParam = SP))
   
     		# invoke gsModel function to fit GBLUP model
    		modelParms <- gsModel(snpsMarker = trainRec$genoNA, pheno  = trainRec$phenoNA)


    		# make new crosses before 100 crosses 50 progeny
    		F1 <- randCross(pop = parents, nCrosses = 100, nProgeny = nProgeny,simParam = SP)

    		# save mean and genetic variance to track progress of breeding program
    		meanGV[year] <- meanG(PYT)
		varGV[year] <- varG_ADG(PYT,pval)

		# heritability estimates
                h2[year] <- modelParms[[2]]
                H2[year] <- modelParms[[4]]

	} # end loop 10-cycle of selection
 
#  put the simulated parameters as a data frame object  
   simParms <- data.frame(simRun=rep(REP,nCycles),
                        year=(1:nCycles),
                        scenario= rep("Narrow adaptation", nCycles),
                        meanGV=meanGV,
                        varGV=varGV,
			varGxE=varGxE,
			h2=h2,
			H2=H2,							
                        accCET=accCET,
                        accPYT=accPYT,
                        accAYT=accAYT,
                        accUYT=accUYT,

                        nParCET=nParCET,
                        nParPYT=nParPYT,
                        nParAYT=nParAYT,
                        nParUYT=nParUYT)

   return(simParms)

} # end gsNarrow function
