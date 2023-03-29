# A function defined to simulate narrow adaptation
gsNarrow  <- function(REP,ME,pvalUYT,locUYT,
	            pvalAYT, locAYT,
		    pvalPYT,locPYT){

  cat("\n Narrow breeding program for mega environment ",ME, "\n")

  cat("\n Simulation run:",REP,"of", nSimRun,"replication", "\n")

  # load the burnin phase for each replicate of this scenario
  load(paste0("./data/burnin_",REP,".rda"))

  if(ME == "ME1"){
        # CET <- splitCET$ME1  # extract the subset of CET for ME1 from TP
        trainRec$phenoNA <- trainRec$phenoME1
        trainRec$genoNA <- trainRec$genoME1
  } else if(ME == "ME2"){

        #CET <- splitCET$ME2  # extract the subset of CET for ME2 from TP
        trainRec$phenoNA <- trainRec$phenoME2
        trainRec$genoNA <- trainRec$genoME2
  }


  # Test the new scenario for 10-cycle of selection after burnin phase


# Generate or replace presample p-value for each cycle year
# in the burn-in and future evaluation steps

allLocPvals <- numeric(9)

pvalRange <- as.data.frame(matrix(c(0.1, 0.5, 0.15, 0.55, 0.2, 0.62, 0.2,0.68,
                      0.3, 0.7, 0.4, 0.8, 0.4, 0.8, 0.45, 0.85, 0.5, 0.9),
                        ncol=2,byrow = TRUE))

for (loc in 1:9){
       allLocPvals[loc] <- runif(1, pvalRange[loc, 1], pvalRange[loc, 2])
}

# invoke gsModel function to fit GBLUP model
gebv <- gsModel(snpsMarker = trainRec$genoNA, datName = trainRec$phenoNA)

# genomic prediction of breeding population at each evaluation stage
    CET@ebv <- as.matrix(gebv[CET@id])
    PYT@ebv <- as.matrix(gebv[PYT@id])
    AYT@ebv <- as.matrix(gebv[AYT@id])
    UYT@ebv <- as.matrix(gebv[UYT@id])


for(year in (burninYears+1):nCycles){ # begin loop of 10-cycle selection
    
    cat("Replication", REP, "Advancing breeding with narrow-adaptation program",ME,
        "year:",year,"of", nCycles, "\n")

    # select variety for release
    variety <- selectInd(pop = UYT,nInd = nVarietySel,trait = 1,use = "pheno",simParam = SP)

  
    # Uniform Yield Trial (UYT)
    UYT <- selectInd(pop = AYT, nInd = nUYT, trait = 1, use = "ebv", simParam = SP)

    # Invoke the function to phenotype selected UYT clones in 4 locations
    UYTrec <- gxeSim(pvalVec = pvalUYT, pop=UYT,
                  varE = errVarUYT, nreps = repUYT,locID = locUYT,year=year)
    
    UYT <- UYTrec[[1]] # extract out UYT pop from the list output
 
  
  #  Advance Yield Trial (AYT)

    AYT <- selectInd(pop = PYT, nInd = nAYT,trait = 1, use = "ebv", simParam = SP)

    # Invoke the function to phenotype selected AYT clones in 2 locations
    AYTrec <- gxeSim(pvalVec = pvalAYT,pop = AYT,
                  varE = errVarAYT, nreps = repAYT,
                  locID = locAYT,year=year)
    
    AYT <- AYTrec[[1]]


    # Preliminary Yield Trial (PYT)    

    PYT <- selectInd(pop = CET, nInd = nPYT, trait = 1, use = "ebv", simParam=SP)   
    PYTrec <- gxeSim(pvalVec = pvalPYT, pop = PYT,
                     varE = errVarPYT,nreps = repPYT,locID = locPYT,year=year)
    PYT <- PYTrec[[1]]
    

    # Clonal Evaluation Trial (CET) - Implement GS
    # Select individuals from the SDN stage based on phenotype performance

    CET <- selectWithinFam(pop = SDN, nInd = famSize, trait = 1, use = "pheno", simParam = SP)
    CET <- selectInd(pop = CET, nInd = nCET, trait = 1, use = "pheno", simParam = SP)

    CETrec <- gxeSim(pvalVec = allLocPvals[5], pop = CET,
                     varE=errVarCET,nreps=repCET,locID = 5,year=year)
    
    CET <- CETrec[[1]] # CET population extract from the list 


    # seedling nursery
    SDNrec  <- gxeSim(pvalVec=allLocPvals[5], pop = F1,
                   varE = errVarSDN,nreps = repSDN, locID = 5,year=year)

    SDN <- SDNrec[[1]]

    # Update training population
    # by retaining the last 2 year of training data

    # Updates the phenotypic data
    temp <- rbind(CETrec[[2]],PYTrec[[2]],AYTrec[[2]],UYTrec[[2]])
    trainRec$phenoNA <- rbind(trainRec$phenoNA[-(1:(nInd(CET)*1+nInd(PYT)*1+nInd(AYT)*2+nInd(UYT)*4)),],temp)

 
    # update genotyping data
    trainRec$genoNA <- rbind(trainRec$genoNA[-(1:nInd(c(CET, PYT,AYT,UYT))),],
                             pullSnpGeno(pop=c(CET, PYT, AYT, UYT),
                                         simParam = SP))


    # invoke gsModel function to fit GBLUP model
    gebv <- gsModel(snpsMarker = trainRec$genoNA, datName = trainRec$phenoNA)

    # genomic prediction of breeding population at each evaluation stage
    CET@ebv <- as.matrix(gebv[CET@id])
    PYT@ebv <- as.matrix(gebv[PYT@id])
    AYT@ebv <- as.matrix(gebv[AYT@id])
    UYT@ebv <- as.matrix(gebv[UYT@id])

    if (any(is.na(CET@ebv))){
      fileName <- paste0("MissingEBV_CET_rep", REP, "_year", year, ".rds")
      saveRDS(mget(ls()), fileName)
      CET <- CET[!is.na(CET@ebv)] # remove missing value
   }

    # Recycle  new parents in the current year based on gebv before advancing materials
    
    tryCatch(
    { # catch error from selecting from parent pools
    parSelCand <- c(UYT,AYT,PYT,CET)
    parSelCand@ebv <- as.matrix(gebv[parSelCand@id])

    if (any(is.na(parSelCand@ebv))){
      fileName <- paste0("Missing_parentEBV_rep", REP, "_year", year, ".rds")
      saveRDS(mget(ls()), fileName)
    }
    parSelCand <- parSelCand[!is.na(parSelCand@ebv)] # remove missing value
    parents <- selectInd(pop = parSelCand, nInd = nParents,trait = 1, use = "ebv",simParam = SP)

    }, # close brace tryCatch()

    error=function(e) { # begin error
        message('An error occurred while selecting from parental candidate pool')
        print(e)
   } # end error
  ) # close tryCatch() for parent


    # number of individuals by stages in the parental candidate
    nParCET[year] <- sum(parents@id %in% CET@id)
    nParPYT[year] <- sum(parents@id %in% PYT@id)
    nParAYT[year] <- sum(parents@id %in% AYT@id)
    nParUYT[year] <- sum(parents@id %in% UYT@id)


    # make new crosses before 100 crosses 50 progeny
    F1 <- randCross(pop = parents, nCrosses = nCrosses, nProgeny = nProgeny,simParam = SP)

    # save mean and genetic variance to track progress of breeding program
    meanGV[year] <- meanG(PYT)
    varGV[year] <- varG(PYT)

    #eGain[year] = (nPYT*cor(gv(PYT), pheno(PYT))*sqrt(varG(PYT)))/4

    H2[year] <- varModel(METdat = UYTrec[[2]]) # estimate H2 at UYT

    # Accuracies of selection from each breeding stage
    accUYT[year] <-  cor(gv(UYT), ebv(UYT))
    accAYT[year] <-  cor(gv(AYT), ebv(AYT))
    accPYT[year] <-  cor(gv(PYT), ebv(PYT))
    accCET[year] <-  cor(gv(CET), ebv(CET))    

  } # end loop 10-cycle of selection
 
#  put the simulated parameters as a data frame object  
   simParms <- data.frame(simRun=rep(REP,nCycles),
                         year=(1:nCycles),
                         scenario= rep("Narrow adaptation", nCycles),
                         meanGV=meanGV,
                         varGV=varGV,
			 varGxE=varGxE,
			 #eGain=eGain,

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
