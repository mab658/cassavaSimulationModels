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



# Generate or replace presample p-value for each cycle year
# in the burn-in and future evaluation steps

allLocPvals <- numeric(9)
pvalRange <- as.data.frame(matrix(c(0.1, 0.5, 0.15, 0.55, 0.2, 0.62, 0.2,0.68,
                           0.3, 0.7, 0.4, 0.8, 0.4, 0.8, 0.45, 0.85, 0.5, 0.9),
                            ncol=2,byrow = TRUE))

for (loc in 1:9){
       allLocPvals[loc] <- runif(1, pvalRange[loc, 1], pvalRange[loc, 2])
}

  # train the GS model for genomic prediction
  gebv <- gsModel(snpsMarker=trainRec$genoBA,datName=trainRec$phenoBA)

  # genomic prediction of breeding population at each evaluation stage
  CET@ebv <- as.matrix(gebv[CET@id])
  PYT@ebv <- as.matrix(gebv[PYT@id])
  AYT@ebv <- as.matrix(gebv[AYT@id])
  UYT@ebv <- as.matrix(gebv[UYT@id])

  for(year in (burninYears+1):nCycles){ # beginning of Genomic Selection
  	cat("Replication", REP, "Advancing breeding with broad-adaptation program
      	strategy year:",year,"of", nCycles, "\n")
  

  # Selecting variety for release from late stage of evaluation
  variety <- selectInd(pop = UYT, nInd = nVarietySel,trait = 1, use = "ebv",simParam = SP)

  UYT <- selectInd(pop = AYT, nInd = nUYT, trait = 1, use = "ebv", simParam = SP)

  # Invoke the function to phenotype selected UYT clones in 8 locations
  
  
  UYTrec <- gxeSim(pvalVec=allLocPvals[-5], pop = UYT,
                varE = errVarUYT, nreps = repUYT, locID=c(1:4,6:9),year=year)
  UYT <- UYTrec[[1]]  



  # Advance Yield Trial (AYT) - use EBV to select the best AYT lines from PYT

  AYT <- selectInd(pop = PYT, nInd = nAYT,trait = 1,use = "ebv",simParam = SP)
 
  # Invoke the function to phenotype selected AYT clones in 4 locations
  AYTrec <- gxeSim(pvalVec=allLocPvals[c(3,4,6,7)],pop = AYT,
                varE = errVarAYT, nreps = repAYT,  locID = c(3,4,6,7),year=year)
  AYT <- AYTrec[[1]]


  # Preliminary  Yield  Trial (PYT) - use EBV to select the best CET to  PYT lines

  PYT <- selectInd(pop = CET, nInd = nPYT, use = "ebv", simParam = SP)

  # Invoke the function to phenotype selected PYT clones in 2 locations
  PYTrec <- gxeSim(pvalVec=allLocPvals[c(3,7)],pop = PYT,
                varE = errVarPYT, nreps = repPYT, locID=c(3,7),year=year)
  PYT <- PYTrec[[1]]  

  
  # Clonal Evaluation Trial (CET) - Implement GS

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
  gebv <- gsModel(snpsMarker=trainRec$genoBA,datName=trainRec$phenoBA)

  # genomic prediction of breeding population at each evaluation stage
  CET@ebv <- as.matrix(gebv[CET@id])
  PYT@ebv <- as.matrix(gebv[PYT@id])
  AYT@ebv <- as.matrix(gebv[AYT@id])
  UYT@ebv <- as.matrix(gebv[UYT@id])

  if (any(is.na(CET@ebv))){
      fileName <- paste0("MissingEBV_CET_rep", REP, "_year", year, ".rds")
      saveRDS(mget(ls()), fileName)
  }

  CET <- CET[!is.na(CET@ebv)] # remove missing value


  # recycle new parent by selecting base on ebv from UYT, AYT, PYT and CET

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


  # number of individuals by  trial stages in the parental candidate
  nParCET[year] <- sum(parents@id %in% CET@id)
  nParPYT[year] <- sum(parents@id %in% PYT@id)
  nParAYT[year] <- sum(parents@id %in% AYT@id)
  nParUYT[year] <- sum(parents@id %in% UYT@id)

  # crossing block
  F1 <- randCross(pop = parents, nCrosses = nCrosses, nProgeny = nProgeny, simParam = SP)


  # to track progress of breeding program at PYT stage
  # Save mean and variance of genetic values
  meanGV[year] <- meanG(PYT)
  varGV[year] <- varG(PYT)

  #eGain[year] = (nPYT*cor(gv(PYT), pheno(PYT))*sqrt(varG(PYT)))/4

 # estimate heritability from MET data from UYT
 H2[year] <- varModel(METdat = UYTrec[[2]]) # estimate H2 at UYT

  # Accuracies of selection from each breeding stage
  accUYT[year] <-  cor(gv(UYT), ebv(UYT))
  accAYT[year] <-  cor(gv(AYT), ebv(AYT))
  accPYT[year] <-  cor(gv(PYT), ebv(PYT))
  accCET[year] <-  cor(gv(CET), ebv(CET))

} # end loop

# put the simulated parameters as a data frame object
  simParms <- data.frame(simRun=rep(REP,nCycles),
                     year=(1:nCycles),
                     scenario= rep("Broad adaptation", nCycles),
                     meanGV=meanGV,
                     varGV=varGV,
                     varGxE=varGxE,
		    # eGain = eGain,
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
