# This script implemented conventional breeding scheme bssed on phenotypic selection
baseMod <- function(REP){

	cat("\n Modeling scenario 1 conventional breeding scheme","\n")

   	cat("\n Simulation run:",REP,"of", nSimRun,"replication", "\n")

   	# load the state of simulation at last year of burnin phase
   	load(paste0("./data/burnin_",REP,".rda"))

	for(year in (burninYears+1):nCycles){
 		# Desired number of years of a future breeding program
  		cat("\n Advance breeding with conventional breeding strategy year:",year,"of", burninYears+futureYears,"\n")

		# pool the data from PYT, AYT, and UYT for heritability estimate

		if(year == (burninYears+1)){

                        # phenotyping indvs from PYT, AYT, and UYT TP
                        gxeDat  <- rbind(PYTrec[[2]], AYTrec[[2]],UYTrec[[2]])
                } else if(year > (burninYears+1) && year <= nCycles){

                        # update TP
                        gxeDat <- rbind(gxeDat,PYTrec[[2]],AYTrec[[2]],UYTrec[[2]])

                } # end else if

	
		#  recycle parents each year for crossing block
                parents <- selectInd(pop=c(UYT,AYT),nInd=nParents, use="pheno",simParam=SP)

                # number of trials by stages in the parental candidate
		nParCET[year] <- NA
                nParPYT[year] <- NA
                nParAYT[year] <- sum(parents@id %in% AYT@id)
                nParUYT[year] <- sum(parents@id %in% UYT@id)

  		# Assign p-value to each location in a cycle year
  		for (loc in 1:9){
       			allLocPvals[loc] <- pnorm(runif(1, eCovRange[loc,1], eCovRange[loc,2]))
  		}

		# select variety for release
                accUYT[year] <-  cor(gv(UYT), pheno(UYT)) 
  		variety <- selectInd(pop = UYT,nInd = nVarietySel,trait = 1, use = "pheno",simParam = SP)
  
		# Uniform Yield Trial (UYT)
                accAYT[year] <-  cor(gv(AYT), pheno(AYT))
  		UYT <- selectInd(pop = AYT, nInd = nUYT,trait = 1, use = "pheno",simParam = SP)
  

  		# Invoke the function to phenotype selected UYT clones in 10 locations

  		UYTrec <- gxeSim(pvalVec=allLocPvals[-5], pop = UYT,
                	varE = errVarUYT, nreps = repUYT,locID=c(1:4,6:9),year=year)
  		UYT <- UYTrec[[1]]


  		# Advance Yield Trial (AYT)
		accPYT[year] <-  cor(gv(PYT), pheno(PYT))
  		AYT <-  selectInd(pop = PYT, nInd = nAYT, trait = 1,use = "pheno",simParam = SP)

  		# Invoke the function to phenotype selected AYT clones in 4 locations
  		AYTrec <- gxeSim(pvalVec=allLocPvals[c(3,4,6,7)],pop = AYT,
                	varE = errVarAYT, nreps = repAYT, locID = c(3,4,6,7),year=year)
  		AYT <- AYTrec[[1]]


  		# Preliminary  Yield  Trial (PYT)
		accCET[year] <-  cor(gv(CET), pheno(CET))
  		PYT <- selectInd(pop = CET, nInd=nPYT,trait = 1, use = "pheno",simParam = SP)

  		# Invoke the function to phenotype selected PYT clones in 2 locations
  		PYTrec <- gxeSim(pvalVec=allLocPvals[c(3,7)],pop = PYT,
                	varE = errVarPYT, nreps = repPYT,locID=c(3,7),year=year)
  		PYT <- PYTrec[[1]]

		# Clonal Evaluation Trial (CET)
		accSDN[year] <-  cor(gv(SDN), pheno(SDN))
  		CET <- selectWithinFam(pop = SDN, nInd = famSize, trait = 1, use = "pheno", simParam=SP)

  		#Invoke the function to phenotype selected SDNclones in 1 locations
  		CETrec <- gxeSim(pvalVec=allLocPvals[5], pop = CET,
                   	varE = errVarCET,nreps = repCET,locID = 5,year=year)
  		CET <- CETrec[[1]]


  		# Seedling Nursery (SDN)
		SDN <- F1
   		SDNrec  <- gxeSim(pvalVec=allLocPvals[5], pop = SDN,
                	varE = errVarSDN,nreps = repSDN, locID = 5,year=year)

  		SDN <- SDNrec[[1]]


  		# Crossing block - Make crosses to generate new population
  		F1 <- randCross(pop = parents, nCrosses = nCrosses,
                  	nProgeny = nProgeny, simParam = SP)

  
  		# save mean & variance of genetic values to track progress of breeding programs
		meanGV[year] <- meanG(CET)
                varGV[year] <- varG(CET)

		h2[year] <- NA # narrow-sense heritability not estimable
                H2[year] <- varModel(METdat = gxeDat) # broad-sense heritability  
	} # end loop

	# put the simulated parameters as a data frame object
	simParms <- data.frame(simRun=rep(REP,nCycles),
 		year=(1:nCycles),
	       	scenario= rep("Conv", nCycles),
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

		nParCET=nParCET,
               	nParPYT=nParPYT,
               	nParAYT=nParAYT,
               	nParUYT=nParUYT)
	return(simParms)
} # end function
