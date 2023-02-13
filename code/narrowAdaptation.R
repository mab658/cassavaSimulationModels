# A function defined to simulate narrow adaptation
gsNarrow  <- function(REP,ME,pval1UYT,pval2UYT,
                    pval1AYT,pval2AYT,
                    pvalPYT){

  cat("\n Narrow breeding program for mega environment ",ME, "\n")

  cat("\n Simulation run:",REP,"of", nSimRun,"replication", "\n")

  # load the burnin phase for each replicate of this scenario
  load(paste0("./data/burnin_",REP,".rda"))

  if(ME == "ME1"){
        CET <- splitCET$ME1  # extract the subset of CET for ME1 from TP
        trainRec$phenoNA <- trainRec$phenoME1
        trainRec$genoNA <- trainRec$genoME1
  } else if(ME == "ME2"){

        CET <- splitCET$ME2  # extract the subset of CET for ME2 from TP
        trainRec$phenoNA <- trainRec$phenoME2
        trainRec$genoNA <- trainRec$genoME2
  }

  # Test the new scenario for 10-cycle of selection after burnin phase

  for(year in (burninYears+1):nCycles){ # begin loop of 10-cycle selection
    
   cat("Advancing breeding with Narrow-adaptation program year:",year,"of", nCycles, "\n")


   # invoke gsModel function to fit GBLUP model
   gebv <- gsModel(snpsMarker = trainRec$genoNA, datName = trainRec$phenoNA)


   # Recycle  new parents in the current year based on gebv before advancing materials

    tryCatch(
  {
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


    # number of individuals by stages in the parental candidate
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


    # select variety for release
    variety <- selectInd(pop = UYT,nInd = nVarietySel,trait = 1,use = "pheno",simParam = SP)
  

    # Uniform Yield Trial (UYT)

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

    # Invoke the function to phenotype selected UYT clones in 4 locations
    UYTrec <- gxeSim(pval1 = pval1UYT,pval2 = pval2UYT, pop=UYT,
                  varE = errVarUYT, nreps = repUYT,nLocs = 4)
    
    UYT <- UYTrec[[1]] # extract out UYT pop from the list output
 


  
    #  Advance Yield Trial (AYT)

    tryCatch(
    { # catch error at AYT
        AYT <- selectInd(pop = PYT, nInd = nAYT,trait = 1, use = "ebv", simParam = SP)

    }, # close brace tryCatch()

     error=function(e) { # begin error
        message('An error occurred while selecting from PYT to AYT')
        print(e)
        saveRDS(mget(ls()),file="missingSelectionPYT_to_AYT.rds")

} # end error
) # close tryCatch() for AYT

    # Invoke the function to phenotype selected AYT clones in 2 locations
    AYTrec <- gxeSim(pval1 = pval1AYT, pval2 = pval2AYT,pop = AYT,
                  varE = errVarAYT, nreps = repAYT,
                  nLocs = 2)
    
    AYT <- AYTrec[[1]]



    # Preliminary Yield Trial (PYT)    

    tryCatch(
    { # catch error at PYT
      PYT <- selectInd(pop = CET, nInd = nPYT, trait = 1, use = "ebv", simParam=SP)   
    }, # close brace tryCatch()

    error=function(e) { # begin error
        message('An error occurred while selecting from CET to PYT')
	saveRDS(mget(ls()),file=paste0("missingSelectionCET_to_PYT",REP,"_",year,".rds"))
        print(e)
   } # end error

   ) # close tryCatch() for PYT

    PYTrec <- gxeSim(pval1 = pvalPYT, pval2 = pvalPYT, pop = PYT,
                     varE = errVarPYT,nreps = repPYT,nLocs = 1)
    PYT <- PYTrec[[1]]
    

    # Clonal Evaluation Trial (CET) - Implement GS
    # Select individuals from the SDN stage based on phenotype performance

    tryCatch(
    { # catch error at CET
      CET <- selectWithinFam(pop = SDN, nInd = 10, trait = 1, use = "pheno", simParam = SP)
      CET <- selectInd(pop = CET, nInd = 500, trait = 1, use = "pheno", simParam = SP)
    }, # close brace tryCatch()

    error=function(e) { # begin error
        message('An error occurred while selecting from SDN to CET')
	saveRDS(mget(ls()),file="missingSelectionSDN_to_CET.rds")
        print(e)

   } # end error
   ) # close tryCatch() for CET

    CETrec <- gxeSim(pval1 = 0.5, pval2 = 0.5, pop = CET,
                     varE=errVarCET,nreps=repCET,nLocs = 1)
    
    CET <- CETrec[[1]] # CET population extract from the list 


    # seedling nursery
    SDN <- setPheno(pop = F1, varE = errVarSDN, reps = repSDN, p = 0.5, simParam = SP)

    # make new crosses
    F1 <- randCross(pop = parents, nCrosses = 100, nProgeny = 50,simParam = SP)


    # save mean and genetic variance to track progress of breeding program
    meanGV[year] <- meanG(PYT)
    varGV[year] <- varG(PYT)


    # Update training population
    # by retaining the last 2 year of training data

     # Updates the phenotypic data
    temp <- rbind(CETrec[[2]],PYTrec[[2]],AYTrec[[2]],UYTrec[[2]])
    trainRec$phenoNA <- rbind(trainRec$phenoNA[-(1:nInd(c(CET, PYT,AYT,UYT))),],temp)


    # update genotyping data
    trainRec$genoNA <- rbind(trainRec$genoNA[-(1:nInd(c(CET, PYT,AYT,UYT))),],
                             pullSnpGeno(pop=c(CET, PYT, AYT, UYT),
                                         simParam = SP))
    
    # remove duplicates individuals from genotyping
    trainRec$genoNA <- trainRec$genoNA[!duplicated(rownames(trainRec$genoNA)),]

  } # end loop 10-cycle of selection
 
 
#  put the simulated parameters as a data frame object  
   simParms <- data.frame(simRun=rep(REP,nCycles),
                         year=(1:nCycles),
                         scenario= rep("Narrow adaptation", nCycles),
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

} # end gsNarrow function
