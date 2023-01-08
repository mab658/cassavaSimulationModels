# Scenario 2
# This script implemented broad adaptation breeding programs
# with GS to estimate BV in CET based on genomic data to boost 
# the selection accuracy at this early stage where phenotype 
# information is limited
# GS is used to advance from CET, PYT, AYT, and UYT

# Train genomic selection (GS) model and safe the model fit

# invoke gsModel function to fit GBLUP model
gebv <- gsModel(snpsMarker=snps,datName=phenoDat)

# Estimate breeding values in  CET, PYT, AYT, and UYT for being genotyped using GBLUP model
  CET@ebv <- as.matrix(gebv[rownames(gebv) %in% CET@id])
  PYT@ebv <- as.matrix(gebv[rownames(gebv) %in% PYT@id])
  AYT@ebv <- as.matrix(gebv[rownames(gebv) %in% AYT@id])
  UYT@ebv <- as.matrix(gebv[rownames(gebv) %in% UYT@id])
  
for(year in (burninYears+1):nCycles){ # beginning of Genomic Selection
  cat("Advancing breeding with broad-adaptation program
      strategy year:",year,"of", nCycles, "\n")

  # Selecting variety for release
  variety <- selectInd(UYT,nInd=nVarietySel, use="pheno",simParam=SP)

 
  # Uniform yield trial  (UYT)
  UYT <- selectInd(pop=AYT,  nInd=nUYT, use="ebv", simParam=SP)
 
  # Invoke the function to phenotype selected UYT clones in 8 locations
  UYTrec <- gxeSim(pval1=0.1, pval2=0.9, pop=UYT,
                varE=errVarUYT, nreps=repUYT,nLocs=8)
  UYT <- UYTrec[[1]]
  
  # estimate breeding values for the evaluated UYT lines and evaluated selection accuracy
  UYT@ebv <- as.matrix(gebv[rownames(gebv) %in% UYT@id])
  accUYT[year] <-  cor(gv(UYT), ebv(UYT))



  # Advance Yield Trial (AYT) - use EBV to select the best AYT lines from PYT
  AYT <- selectInd(pop=PYT, nInd=nAYT,use="ebv",simParam=SP)
 
  # Invoke the function to phenotype selected AYT clones in 4 locations
  AYTrec <- gxeSim(pval1=0.3, pval2=0.7,pop=AYT,
                varE=errVarAYT, nreps=repAYT,nLocs=4)
  AYT <- AYTrec[[1]]

  # estimate breeding values for the evaluated AYT lines and evaluated selection accuracy
  AYT@ebv <- as.matrix(gebv[rownames(gebv) %in% AYT@id])
  accAYT[year] <-  cor(gv(AYT), ebv(AYT))


  # Preliminary  Yield  Trial (PYT) - use EBV to select the best PYT lines from CET
  PYT <- selectInd(pop=CET, nInd=nPYT, use="ebv",simParam=SP)
 
  # Invoke the function to phenotype selected PYT clones in 2 locations
  PYTrec <- gxeSim(pval1=0.3, pval2=0.7, pop=PYT,
                varE=errVarPYT, nreps=repPYT, nLocs=2)
  PYT <- PYTrec[[1]]  
  # estimate breeding values for the evaluated PYT  lines 
  PYT@ebv <- as.matrix(gebv[rownames(gebv) %in% PYT@id])
  accPYT[year] <-  cor(gv(PYT), ebv(PYT)) # evaluate accuracy of selection based on GEBV
  

  # Clonal Evaluation Trial (CET) - Implement GS
  # Select individuals from the SDN stage based on phenotype performance
  cat("\n CET corr in loop:",REP,"of", nSimRun,"replication", "\n")
  accCET[year] <-  cor(gv(CET), ebv(CET)) # evaluate accuracy of selection base>
  
  CET <- selectWithinFam(pop=SDN, nInd=famSize, use="pheno",simParam=SP)
  CET <- selectInd(pop=CET, nInd=nCET, use="pheno", simParam=SP)
  #Invoke the function to phenotype selected PYT clones in 2 locations
  CETrec <- gxeSim(pval1=0.5, pval2=0.5, pop=CET,
                   varE=errVarCET,nreps=repCET,nLocs=1)  
  CET <- CETrec[[1]]

  # Seedling Nursery
  SDN <- setPheno(pop=F1, varE=errVarSDN, reps=repSDN, fixEff=year,
                  p=0.5, simParam=SP)



  # Save mean and variance of genetic values
  # to track progress of breeding program
  meanGV[year] <- meanG(PYT)
  varGV[year] <- varG(PYT)

  # Update training population and refit GBLUP model
  # by retaining the last 2 year of training data
  
   # Updates the phenotypic data
  temp <- rbind(CETrec[[2]],PYTrec[[2]],
                AYTrec[[2]],UYTrec[[2]])
  trainRec$phenoBA <- rbind(trainRec$phenoBA[-(1:nInd(c(CET, PYT,AYT,UYT))),],
                            temp)


  # update genotyping data
  trainRec$genoBA <- rbind(trainRec$genoBA[-(1:nInd(c(CET, PYT,AYT,UYT))),],
                           pullSnpGeno(pop=c(CET, PYT, AYT, UYT),
                                       simParam=SP))
  # remove duplicates individuals from genotyping
  trainRec$genoBA <- trainRec$genoBA[!duplicated(rownames(trainRec$genoBA)),]

  phenoDat <- trainRec$phenoBA
  snps <- trainRec$genoBA  

  # Train genomic selection (GS) model and safe the model fit
  gebv <- gsModel(snpsMarker=snps,datName=phenoDat) # invoke gsModel function to fit GBLUP model

  # Estimate breeding values in  CET, PYT, AYT, and UYT for being genotyped using GBLUP model
  CET@ebv <- as.matrix(gebv[rownames(gebv) %in% CET@id])
  PYT@ebv <- as.matrix(gebv[rownames(gebv) %in% PYT@id])
  AYT@ebv <- as.matrix(gebv[rownames(gebv) %in% AYT@id])
  UYT@ebv <- as.matrix(gebv[rownames(gebv) %in% UYT@id])

  # recycle new parent
  parSelCand <- c(UYT,AYT,PYT,CET)
  parSelCand@ebv <- as.matrix(gebv[rownames(gebv) %in% parSelCand@id])
  parents <- selectInd(pop=parSelCand, nInd=nInd(parents), use="ebv",simParam=SP)

  # number of trials by stages in the parental candidate
  nParCET[year] <- sum(parents@id %in% CET@id)
  nParPYT[year] <- sum(parents@id %in% PYT@id)
  nParAYT[year] <- sum(parents@id %in% AYT@id)
  nParUYT[year] <- sum(parents@id %in% UYT@id)

  # crossing block

  F1 <- randCross(pop=parents, nCrosses=nCrosses, nProgeny=nProgeny, simParam=SP)

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
