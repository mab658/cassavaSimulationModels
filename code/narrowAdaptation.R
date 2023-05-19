# A function defined to simulate narrow adaptation
# Test the new breeding scenario for 10-cycle of selection after burnin phase
gsNarrow  <- function(REP,ME,pvalUYT,locUYT,
	            pvalAYT, locAYT,
		    pvalPYT,locPYT){

	cat("\n Narrow breeding program for mega environment ",ME, "\n")

  	cat("\n Simulation run:",REP,"of", nSimRun,"replication", "\n")

  	# load the burnin phase for each replicate of this scenario
  	load(paste0("./data/burnin_",REP,".rda"))


	for(year in (burninYears+1):nCycles){ # begin loop of 10-cycle selection

		cat("Replication", REP, "Advancing breeding with narrow-adaptation program",ME,
                        "year:",year,"of", nCycles, "\n")


		# Constitute the training population (TP) using PYT, AYT, and UYT
                # with a sliding-window process of 5 year.
                # Accumulation of the records starts a year after  burn-in year
                # the sliding-window process removes the oldest data
                # The number of years retained is 5.

		if(ME == "ME1"){
			pval=0.64
			if(year == (burninYears+1)){
                        	# phenotyping and genotyping indvs in TP set(PYT, AYT, UYT)
 		               	trainRec$phenoNA <- rbind(PYTrec[[2]][PYTrec[[2]]$loc=="L7",],
                       				AYTrec[[2]][AYTrec[[2]]$loc %in% c("L6","L7"),],
                       				UYTrec[[2]][UYTrec[[2]]$loc %in% c("L6","L7","L8","L9"),])

                		trainRec$genoNA <- pullSnpGeno(pop = c(PYT, AYT, UYT), simParam = SP)

                	} else if(year > (burninYears+1) && year <= nCycles){
                        	# update Training population
				trainRec$phenoNA <- rbind(trainRec$phenoNA,
                       		PYTrec[[2]][PYTrec[[2]]$loc=="L7",],
                       		AYTrec[[2]][AYTrec[[2]]$loc %in% c("L6","L7"),],
                       		UYTrec[[2]][UYTrec[[2]]$loc %in% c("L6","L7","L8","L9"),])

	                	trainRec$genoNA <- rbind(trainRec$genoNA,pullSnpGeno(pop=c(PYT, AYT, UYT),simParam=SP))
                	} # end elseif


		} else if(ME == "ME2"){ # Mega-environmet 2
			pval=0.36
			if(year == (burninYears+1)){
                                # phenotyping and genotyping indvs in TP set (PYT,AYT, UYT)
                                trainRec$phenoNA <- rbind(PYTrec[[2]][PYTrec[[2]]$loc=="L3",],
                                                AYTrec[[2]][AYTrec[[2]]$loc %in% c("L3","L4"),],
                                                UYTrec[[2]][UYTrec[[2]]$loc %in% c("L1","L2","L3","L4"),])

                                trainRec$genoNA <- pullSnpGeno(pop = c(PYT, AYT, UYT), simParam = SP)

			} else if(year > (burninYears+1) && year <= nCycles){
                       
                        	# update the training population
                        	trainRec$phenoNA <- rbind(trainRec$phenoNA,
                        	PYTrec[[2]][PYTrec[[2]]$loc=="L3",],
                        	AYTrec[[2]][AYTrec[[2]]$loc %in% c("L3","L4"),],
                        	UYTrec[[2]][UYTrec[[2]]$loc %in% c("L1","L2","L3","L4"),])

                        	trainRec$genoNA <- rbind(trainRec$genoNA,pullSnpGeno(pop=c(PYT, AYT, UYT),simParam=SP))

			} # end elseif

        	} # end if Mega-environment 2

		# invoke gsModel function to fit GBLUP model
        	modelParms <- gsModel(snpsMarker = trainRec$genoNA, pheno  = trainRec$phenoNA)

		# Select  new parents in the current year based on gebv for recycling
                # Use GP model to assign gebv to individuals in the breeding population

                PYT@ebv <- as.matrix(modelParms[[1]][PYT@id])
                AYT@ebv <- as.matrix(modelParms[[1]][AYT@id])
                UYT@ebv <- as.matrix(modelParms[[1]][UYT@id])

		if (any(is.na(PYT@ebv))){
                        #fileName <- paste0("Missing_PYT_rep", REP, "_year", year, ".rds")
                        #saveRDS(mget(ls()), fileName)
                        PYT <- PYT[!is.na(PYT@ebv)] # remove missing value
                }


                # constitute parental selection candidate and select based on GEBV for next cycle
                parSelCand <- c(UYT,AYT,PYT)

                if (any(is.na(parSelCand@ebv))){
                        fileName <- paste0("Missing_parentEBV_rep", REP, "_year", year, ".rds")
                        saveRDS(mget(ls()), fileName)
                        parSelCand <- parSelCand[!is.na(parSelCand@ebv)] # remove missing value
                }

                parents <- selectInd(pop = parSelCand, nInd = nParents,trait = 1, use = "ebv",simParam = SP)

                # number of individuals constitute parental candidates per trial stage
                nParPYT[year] <- sum(parents@id %in% PYT@id)
                nParAYT[year] <- sum(parents@id %in% AYT@id)
                nParUYT[year] <- sum(parents@id %in% UYT@id)


		# Assign GETGV to the population from GP model as surrogates of selection values
        	#PYT@ebv <- as.matrix(modelParms[[3]][PYT@id])
        	#AYT@ebv <- as.matrix(modelParms[[3]][AYT@id])
        	#UYT@ebv <- as.matrix(modelParms[[3]][UYT@id])


		# assign p-value as environmental covariate to each location
  		for (loc in 1:9){
       			allLocPvals[loc] <- runif(1, pvalRange[loc, 1], pvalRange[loc, 2])
  		}
 
		# Selecting from UYT lines based on GETGV for variety for release and its accuracy
		accUYT[year] <- cor(gvADG(UYT,pval), ebv(UYT))

    		variety <- selectInd(pop = UYT,nInd = nVarietySel,trait = 1,use = "ebv",simParam = SP)
		# Uniform Yield Trial (UYT) - use GETGV to select the best from AYT lines to UYT

        	accAYT[year] <- cor(gvADG(AYT,pval),ebv(AYT))
		UYT <- selectInd(pop = AYT, nInd = nUYT, trait = 1, use = "ebv", simParam = SP)

    		# Invoke the function to phenotype selected UYT clones in 4 locations
    		UYTrec <- gxeSim(pvalVec = pvalUYT, pop = UYT,
                  	varE = errVarUYT, nreps = repUYT,locID = locUYT,year=year)
   
    		UYT <- UYTrec[[1]] # extract out UYT pop from the list output
 
  
		#  Advance Yield Trial (AYT) - use GETGV to select the best PYT lines and its accuracy
        	accPYT[year] <- cor(gvADG(PYT,pval), ebv(PYT))
    		AYT <- selectInd(pop = PYT, nInd = nAYT,trait = 1, use = "ebv", simParam = SP)
    
		# Invoke the function to phenotype selected AYT clones in 2 locations
    		AYTrec <- gxeSim(pvalVec = pvalAYT,pop = AYT,
                  	varE = errVarAYT, nreps = repAYT,
                  	locID = locAYT,year=year)
    
    		AYT <- AYTrec[[1]]


		# Preliminary Yield Trial (PYT) - use GETGV to select the best CET lines to  PYT
        	accCET[year] <- cor(gvADG(CET,pval), pheno(CET))
    		PYT <- selectInd(pop = CET, nInd = nPYT, trait = 1, use = "pheno", simParam=SP)   
    		PYTrec <- gxeSim(pvalVec = pvalPYT, pop = PYT,
                   	varE = errVarPYT,nreps = repPYT,locID = locPYT,year=year)
  

    		# Clonal Evaluation Trial (CET) - Implement GS
    		# Select individuals from the SDN stage based on phenotype performance

                accSDN[year] <- cor(gvADG(SDN,pval), pheno(SDN))
    		CET <- selectWithinFam(pop = SDN, nInd = famSize, trait = 1, use = "pheno", simParam = SP)
    		CETrec <- gxeSim(pvalVec = allLocPvals[5], pop = CET,
                     	varE=errVarCET,nreps=repCET,locID = 5,year=year)
    
    		CET <- CETrec[[1]] # CET population extract from the list 


    		# seedling nursery
		SDN <- F1
    		SDNrec  <- gxeSim(pvalVec=allLocPvals[5], pop = SDN,
                	   varE = errVarSDN,nreps = repSDN, locID = 5,year=year)

    		SDN <- SDNrec[[1]]
		
		# Crossing block - random mating of parents to make 100 crosses with 50 progeny per cross
        	F1 <- randCross(pop = parents, nCrosses = 100, nProgeny = nProgeny,simParam = SP)

    		# save mean and genetic variance to track progress of breeding program
		meanGV[year] <- meanG(PYT)
		varGV[year] <- varG_ADG(PYT,pval)

		# heritability estimates
        	h2[year] <- modelParms[[2]] # Narrow-sense heritability
        	H2[year] <- modelParms[[4]] # broad-sense heritability
    
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

		accSDN=accSDN,
         	accCET=accCET,
         	accPYT=accPYT,
         	accAYT=accAYT,
         	accUYT=accUYT,

         	nParPYT=nParPYT,
 		nParAYT=nParAYT,
         	nParUYT=nParUYT)

   	return(simParms)

} # end gsNarrow function
