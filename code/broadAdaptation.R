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
  

  for(year in (burninYears+1):nCycles){ # beginning of Genomic Selection
  	cat("Advancing breeding with broad-adaptation program
      	strategy year:",year,"of", nCycles, "\n")


  # train the GS model for genomic prediction
  gebv <- gsModel(snpsMarker=trainRec$genoBA,datName=trainRec$phenoBA)

  # recycle new parent by selecting base on ebv from UYT, AYT, PYT and CET
  
  tryCatch(
  { # catch error from selecting from parent pools
    parSelCand <- c(UYT,AYT,PYT,CET)
    parSelCand@ebv <- as.matrix(gebv[parSelCand@id])
    parents <- selectInd(pop = parSelCand, nInd = nParents,trait = 1, use = "ebv",simParam = SP)
  }, # close brace tryCatch()

  error=function(e) { # begin error
        message('An error occurred while selecting from parental candidate pool')
        saveRDS(mget(ls()),file=paste0("missingSelectionfrom_parentPool",REP,"_",year,".rds"))
	print(e)
  } # end error
) # close tryCatch() for parent


  # number of trials by stages in the parental candidate
  nParCET[year] <- sum(parents@id %in% CET@id)
  nParPYT[year] <- sum(parents@id %in% PYT@id)
  nParAYT[year] <- sum(parents@id %in% AYT@id)
  nParUYT[year] <- sum(parents@id %in% UYT@id)

  # genomic prediction of breeding population at each evaluation stage
  CET@ebv <- as.matrix(gebv[CET@id])
  PYT@ebv <- as.matrix(gebv[PYT@id])
  AYT@ebv <- as.matrix(gebv[AYT@id])
  UYT@ebv <- as.matrix(gebv[UYT@id])


  # error tracking for population without ebv assigned
   if(any(is.na(CET@ebv)) | any(is.na(PYT@ebv)) |
        any(is.na(AYT@ebv)) | any(is.na(UYT@ebv))){ # begin if
         print("There are missing values for EBVs")
         saveRDS(mget(ls()),file="missingEBVValues.rds")
    } # end if

  # Accuracies of selection from each breeding stage
  accUYT[year] <-  cor(gv(UYT), ebv(UYT))
  accAYT[year] <-  cor(gv(AYT), ebv(AYT))
  accPYT[year] <-  cor(gv(PYT), ebv(PYT))
  accCET[year] <-  cor(gv(CET), ebv(CET))


  # Selecting variety for release from late stage of evaluation
  variety <- selectInd(pop = UYT,nInd = nVarietySel,trait = 1, use = "pheno",simParam = SP)

  tryCatch(
  { # catch error at UYT  
    UYT <- selectInd(pop = AYT, nInd = nUYT, trait = 1, use = "ebv", simParam = SP)
 
  }, # close brace tryCatch()

  error=function(e) { # begin error
        message('An error occurred while selecting from AYT to UYT')
	saveRDS(mget(ls()),file="missingSelection_AYT_to_UYT.rds")
        print(e)

} # end error
) # close tryCatch() for UYT

  # Invoke the function to phenotype selected UYT clones in 8 locations
  UYTrec <- gxeSim(pval1 = 0.1, pval2 = 0.9, pop = UYT,
                varE = errVarUYT, nreps = repUYT,nLocs = 8)
  UYT <- UYTrec[[1]]


  # Advance Yield Trial (AYT) - use EBV to select the best AYT lines from PYT
  tryCatch(
  { # catch error at AYT

     AYT <- selectInd(pop = PYT, nInd = nAYT,trait = 1,use = "ebv",simParam = SP)
  }, # close brace tryCatch()

  error=function(e) { # begin error
        message('An error occurred while selecting from PYT to AYT')
	saveRDS(mget(ls()),file="missingSelectionPYT_to_AYT.rds")
        print(e)
  } # end error

) # close tryCatch() for AYT 

  # Invoke the function to phenotype selected AYT clones in 4 locations
  AYTrec <- gxeSim(pval1 = 0.3, pval2 = 0.7,pop = AYT,
                varE = errVarAYT, nreps = repAYT,nLocs = 4)
  AYT <- AYTrec[[1]]



  # Preliminary  Yield  Trial (PYT) - use EBV to select the best CET to  PYT lines

  tryCatch(
  { # catch error at PYT

    PYT <- selectInd(pop = CET, nInd = nPYT, use = "ebv", simParam = SP)
    cat("CET EBV values at rep ", REP,"year",year,"\n")
    print(CET@ebv)
   }, # close brace tryCatch()

  error=function(e) { # begin error
        message('An error occurred while selecting from CET to PYT')
	saveRDS(mget(ls()),file=paste0("missingSelectionCET_to_PYT",REP,"_",year,".rds"))
	 print(e)
  } # end error

) # close tryCatch() for PYT

  # Invoke the function to phenotype selected PYT clones in 2 locations
  PYTrec <- gxeSim(pval1 = 0.3, pval2 = 0.7, pop = PYT,
                varE = errVarPYT, nreps = repPYT, nLocs = 2)
  PYT <- PYTrec[[1]]  
  

  # Clonal Evaluation Trial (CET) - Implement GS
  tryCatch(
  { # catch error at CET

   CET <- selectWithinFam(pop = SDN, nInd = famSize, trait = 1, use = "pheno", simParam = SP)
   CET <- selectInd(pop = CET, nInd = nCET, trait = 1, use = "pheno", simParam = SP)
  }, # close brace tryCatch()

  error=function(e) { # begin error
        message('An error occurred while selecting from SDN to CET')
	saveRDS(mget(ls()),file="missingSelectionSDN_to_CET.rds")
        print(e)
} # end error
) # close tryCatch() for CET

  #Invoke the function to phenotype selected PYT clones in 2 locations
  CETrec <- gxeSim(pval1 = 0.5, pval2 = 0.5, pop = CET,
                   varE = errVarCET,nreps = repCET,nLocs = 1)  
  CET <- CETrec[[1]]

  # Seedling Nursery
  SDN <- setPheno(pop = F1, varE = errVarSDN, reps = repSDN, fixEff = year,
                  p = 0.5, simParam = SP)

  # crossing block
  F1 <- randCross(pop = parents, nCrosses = nCrosses, nProgeny = nProgeny, simParam = SP)


  # to track progress of breeding program at PYT stage
  # Save mean and variance of genetic values
  meanGV[year] <- meanG(PYT)
  varGV[year] <- varG(PYT)


  # Update training population
  # by retaining the last 2 year of training data
  
  # Updates the phenotypic data
  temp <- rbind(CETrec[[2]],PYTrec[[2]], AYTrec[[2]], UYTrec[[2]])                   
  trainRec$phenoBA <- rbind(trainRec$phenoBA[-(1:nInd(c(CET, PYT,AYT,UYT))),],
                            temp)
   
  # update genotyping data
  trainRec$genoBA <- rbind(trainRec$genoBA[-(1:nInd(c(CET, PYT,AYT,UYT))),],
  pullSnpGeno(pop = c(CET, PYT, AYT, UYT), simParam = SP))

  # remove duplicates individuals from genotyping data
  trainRec$genoBA <- trainRec$genoBA[!duplicated(rownames(trainRec$genoBA)),]

} # end loop

# put the simulated parameters as a data frame object
  simParms <- data.frame(simRun=rep(REP,nCycles),
                     year=(1:nCycles),
                     scenario= rep("Broad adaptation", nCycles),
                     meanGV=meanGV,
                     varGV=varGV,
                     varGxE=varGE,

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
