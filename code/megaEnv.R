megaEnv <- function(pval1UYT,pval2UYT,
                    pval1AYT,pval2AYT,
                    pvalPYT){


# invoke gsModel function to fit GBLUP model
gebv <- gsModel(snpsMarker=trainRec$genoNA,datName=trainRec$phenoNA)


  # assign ebv to the population from GBLUP model
  # assign ebv to the population from GBLUP model
  CET@ebv <- as.matrix(gebv[CET@id])
  PYT@ebv <- as.matrix(gebv[PYT@id])
  AYT@ebv <- as.matrix(gebv[AYT@id])
  UYT@ebv <- as.matrix(gebv[UYT@id])

  for(year in (burninYears+1):nCycles){ # begining of loop
    
    cat("Advancing breeding with Narrow-adaptation program
       year:",year,"of", nCycles, "\n")


    # Selection accuracy from UYT lines for variety release
    accUYT[year] <-  cor(gv(UYT), ebv(UYT))
    variety <- selectInd(pop=UYT,nInd=nVarietySel,use="pheno",simParam=SP)
  

    # Uniform Yield Trial (UYT)
    # Selection accuracy from AYT lines to UYT
    accAYT[year] <-  cor(gv(AYT), ebv(AYT))
   
    UYT <- selectInd(pop=AYT, nInd=nUYT, use="ebv", simParam=SP)
    # Invoke the function to phenotype selected UYT clones in 4 locations
    UYTrec <- gxeSim(pval1=pval1UYT,pval2=pval2UYT, pop=UYT,
                  varE=errVarUYT, nreps=repUYT,nLocs=4)
    
    UYT <- UYTrec[[1]] # extract out UYT pop from the list output
    
  
    #  Advance Yield Trial (AYT)
    # Selection accuracy from PYT lines  to AYT
    accPYT[year] <-  cor(gv(PYT), ebv(PYT))
    AYT <- selectInd(pop=PYT, nInd=nAYT,use="ebv", simParam=SP)

    # Invoke the function to phenotype selected AYT clones in 2 locations
    AYTrec <- gxeSim(pval1=pval1AYT, pval2=pval2AYT,pop=AYT,
                  varE=errVarAYT, nreps=repAYT,
                  nLocs=2)
    
    AYT <- AYTrec[[1]]
    # estimate breeding values for the evaluated AYT lines and evaluated selection accuracy

   
    # Preliminary  Yield  Trial (PYT)
    # evaluate accuracy of selection on CET based on GEBV
    accCET[year] <-  cor(gv(CET), ebv(CET)) # evaluate accuracy of selection based on GEBV
    
    PYT <- selectInd(pop=CET, nInd=nPYT, use="ebv", simParam=SP)
    PYTrec <- gxeSim(pval1=pvalPYT, pval2=pvalPYT, pop=PYT,
                     varE=errVarPYT,nreps=repPYT,nLocs=1)
    PYT <- PYTrec[[1]]
    

    # Clonal Evaluation Trial (CET) - Implement GS
    # Select individuals from the SDN stage based on phenotype performance
    
    # Select individuals from the SDN stage based on phenotype performance
    CET <- selectWithinFam(pop=SDN, nInd=famSize, use="pheno",simParam=SP)
    CET <- selectInd(pop=CET, nInd=500, use="pheno", simParam=SP)
    CETrec <- gxeSim(pval1=0.5, pval2=0.5, pop=CET,
                     varE=errVarCET,nreps=repCET,nLocs=1)
    
    CET <- CETrec[[1]]

    # seedling nursery
    SDN <- setPheno(pop=F1, varE=errVarSDN, reps=repSDN, p=0.5, simParam=SP)
    

    # save mean and genetic variance to track progress of breeding program
    meanGV[year] <- meanG(PYT)
    varGV[year] <- varG(PYT)

    # Update training population to refit GBLUP model
    # by retaining the last 2 year of training data

     # Updates the phenotypic data
    temp <- rbind(CETrec[[2]],PYTrec[[2]],
                  AYTrec[[2]],UYTrec[[2]])
    trainRec$phenoNA <- rbind(trainRec$phenoNA[-(1:nInd(c(CET, PYT,AYT,UYT))),],
                              temp)

    # update genotyping data
    trainRec$genoNA <- rbind(trainRec$genoNA[-(1:nInd(c(CET, PYT,AYT,UYT))),],
                             pullSnpGeno(pop=c(CET, PYT, AYT, UYT),
                                         simParam=SP))
    
    # remove duplicates individuals from genotyping
    trainRec$genoNA <- trainRec$genoNA[!duplicated(rownames(trainRec$genoNA)),]
   
    
    # Train genomic selection (GS) model and safe the model fit
    gebv <- gsModel(snpsMarker=trainRec$genoNA,datName=trainRec$phenoNA) # invoke gsModel function to fit GBLUP model


    # assign ebv to the population from GBLUP model
    CET@ebv <- as.matrix(gebv[CET@id])
    PYT@ebv <- as.matrix(gebv[PYT@id])
    AYT@ebv <- as.matrix(gebv[AYT@id])
    UYT@ebv <- as.matrix(gebv[UYT@id])

    # Recycle  new parents in the current year
    # based on gebv before advancing materials
    parSelCand <- c(UYT,AYT,PYT,CET)
    # parSelCand@ebv <- as.matrix(gebv[parSelCand@id,])
    parents <- selectInd(pop=parSelCand, nInd=nInd(parents), use="ebv",simParam=SP)

    # number of individuals by stages in the parental candidate
    nParCET[year] <- sum(parents@id %in% CET@id)
    nParPYT[year] <- sum(parents@id %in% PYT@id)
    nParAYT[year] <- sum(parents@id %in% AYT@id)
    nParUYT[year] <- sum(parents@id %in% UYT@id)  

    # make new crosses
    F1 <- randCross(pop=parents, nCrosses=100,nProgeny=50,simParam=SP)

  } # end loop
  
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

} # end function
