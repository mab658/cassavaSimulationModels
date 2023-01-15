# After populating the cassava breeding pipeline
# we implement the burn-in phase by advancing population of individuals
# through the breeding scheme over years where the parents are recycled.
# Burnin phase serves as common starting point for
# different breeding scenarios (programs)
# Advance breeding by year by working backwards through
# pipeline  to avoid copying data

# create a placeholder of training records for genomic selection (GS)

trainRec <- vector(mode="list",length=8)
names(trainRec) <- c("phenoBA","genoBA",
                     "phenoME1","genoME1",
                     "phenoME2","genoME2",
                    "phenoNA","genoNA")


# we track the progress of breeding program
# by initializing variables to save the genetic parameters
simRun <- vector("numeric",length=nCycles)
meanGV <- vector("numeric",length=nCycles)
varGV <- vector("numeric",length=nCycles)

accCET <- accPYT <- accAYT <- accUYT <- vector("numeric",length=nCycles)

# pipeline  to avoid copying data

# track number of trial stages in parent per year
nParCET <- nParPYT <- nParAYT <- nParUYT <- vector("numeric",length=nCycles)

# Run 10 years of burn-in (Cycle years) as a common starting point
# for different scenarios.
# This marks the beginning of a new breeding cycle where new parents are
# formed and taking to crossing block

# Advance breeding program by Working backwards through pipeline
# to avoid copying data

for(year in 1:burninYears){#Change to any number of desired years
  cat("Burnin phase cycle year:",year,"of",burninYears, "\n")

  # Select variety to be released

  variety <- selectInd(pop=UYT,nInd=nVarietySel, use="pheno",simParam=SP)

  # Uniform Yield Trial (UYT)

  UYT <- selectInd(pop=AYT, nInd=nUYT, use="pheno",simParam=SP)
  UYTrec <- gxeSim(pval1=0.1, pval2=0.9, pop=UYT,
                varE=errVarUYT, nreps=repUYT, nLocs=8)
 
  UYT <- UYTrec[[1]]

  # Advance Yield Trial (AYT)
  AYT <- selectInd(pop=PYT, nInd=nAYT, use="pheno",simParam=SP)
  AYTrec <- gxeSim(pval1=0.3, pval2=0.7, pop=AYT,
                varE=errVarAYT, nreps=repAYT,nLocs=4)
  
  AYT <- AYTrec[[1]] # extract out the phenotypic values


  # Preliminary  Yield  Trial (PYT)
  PYT <- selectInd(pop=CET, nInd=nPYT, use="pheno", simParam=SP)
  PYTrec <- gxeSim(pval1=0.3, pval2=0.7, pop=PYT,
                varE=errVarPYT,nreps=repPYT,nLocs=2)
  PYT <- PYTrec[[1]]


  #  Clonal Evaluation Trial (CET)
  CET <- selectWithinFam(pop=SDN, nInd=famSize, use="pheno", simParam=SP)
  CET <- selectInd(pop=CET, nInd=nCET, use="pheno", simParam=SP)
  CETrec <- gxeSim(pval1=0.5, pval2=0.5, pop=CET,
                varE=errVarCET,nreps=repCET,nLocs=1)

  CET <- CETrec[[1]]

  # Seedling Nursery
  SDN <- setPheno(pop=F1, varE=errVarSDN,reps=repSDN,fixEff=year,
                  p=0.5,simParam=SP)

  
  # recycle new parents in each year of burn-in phase
  # by selecting best individuals from joint (UYT and AYT) population
 
  parents <- selectInd(c(UYT,AYT), nInd=nInd(parents), 
                       use="pheno",simParam=SP)

  # crossing block
  F1 <- randCross(pop=parents, nCrosses=nCrosses,
                  nProgeny=nProgeny, simParam=SP)

  
  # report the mean and variance of estimated GV in each year
  meanGV[year] <- meanG(PYT)
  varGV[year] <- varG(PYT)
 
  # Constitute the training population (TP)
  # from 8-year cycle to 10-year of burn-in phase
  # using CET,PYT, AYT, and UYT population to have 3 years of TP

  # NOTE: Training records are collected from CET,PYT,AYT,and UYT with
  # a sliding-window process.
  # Accumulation of the records starts in the burn-in year 'startTrainPop'.
  # Once the burn-in period is over, the sliding-window process removes the oldest
  # records. The number of years retained is set with "limityear" in GlobalParameters.

  if(year == startTrainPop){
    # collect phenotype and genotype data of individuals across stages for TP

    # collect phenotype of individuals across stages for TP in broad-adaptation
    trainRec$phenoBA <- rbind(CETrec[[2]],PYTrec[[2]],
                              AYTrec[[2]],UYTrec[[2]])

   
    trainRec$genoBA <- pullSnpGeno(pop=c(CET, PYT, AYT, UYT), 
                                 simParam=SP)
    
    
    # Split the CET to constitute TP for mega-environments- ME1 and ME2
    # for narrow adaptation program
    splitCET <- splitPop(CET)
    
    # collect genotype & phenotype of inds across stages for ME1
    CETrecME1 <-  CETrec[[2]][CETrec[[2]]$id %in% splitCET$ME1@id,]

    trainRec$phenoME1 <- rbind(CETrecME1,PYTrec[[2]],
                              AYTrec[[2]],UYTrec[[2]])

    # genotyping data at beginning of TP for ME1    

    trainRec$genoME1 <- pullSnpGeno(pop=c(splitCET$ME1, 
                                          PYT, AYT, UYT), 
                                 simParam=SP)
    
       
    # collect phenotype and genotype data of individuals for ME2 
    
    CETrecME2 <-  CETrec[[2]][CETrec[[2]]$id %in% splitCET$ME2@id,]
    trainRec$phenoME2 <- rbind(CETrecME2,PYTrec[[2]],
                               AYTrec[[2]],UYTrec[[2]])

    # genotyping data at beginning of TP for ME2
    trainRec$genoME2 <- pullSnpGeno(pop=c(splitCET$ME2, 
                                          PYT, AYT, UYT), 
                                    simParam=SP)
    
    
  } else if(year > startTrainPop && year <= burninYears){
    # update TP for broad adaptation program
    temp <- rbind(CETrec[[2]],PYTrec[[2]],
                  AYTrec[[2]],UYTrec[[2]])
    
    # update TP from year 9 to  10 of burnin phase

    trainRec$phenoBA <- rbind(trainRec$phenoBA,temp)
    trainRec$genoBA <- rbind(trainRec$genoBA,
                             pullSnpGeno(pop=c(CET, PYT, AYT, UYT),
                                         simParam=SP))
    
    # update phenotype and genotype for ME1 and ME2
    splitCET <- splitPop(CET) # Split the CET to constitute TP for ME1 and ME2

    CETrecME1 <-  CETrec[[2]][CETrec[[2]]$id %in% splitCET$ME1@id,]
    trainRec$phenoME1 <- rbind(trainRec$phenoME1,CETrecME1,PYTrec[[2]],
                               AYTrec[[2]],UYTrec[[2]])


    trainRec$genoME1 <- rbind(trainRec$genoME1,
                              pullSnpGeno(pop=c(splitCET$ME1, PYT, AYT, UYT),
                                         simParam=SP))
    
    # update TP for ME2
    
    CETrecME2 <-  CETrec[[2]][CETrec[[2]]$id %in% splitCET$ME2@id,]

    trainRec$phenoME2 <- rbind(trainRec$phenoME2,CETrecME2,PYTrec[[2]],
                               AYTrec[[2]],UYTrec[[2]])

      
    trainRec$genoME2 <- rbind(trainRec$genoME2,
                              pullSnpGeno(pop=c(splitCET$ME2, PYT, AYT, UYT),
                                          simParam=SP))
    
  } # end of elseif
} # end 10-year burn-in

# remove duplicated individuals from genotyping data
trainRec$genoBA <- trainRec$genoBA[!duplicated(rownames(trainRec$genoBA)),]
trainRec$genoME1 <- trainRec$genoME1[!duplicated(rownames(trainRec$genoME1)),]
trainRec$genoME2 <- trainRec$genoME2[!duplicated(rownames(trainRec$genoME2)),]

# Save the state of simulation at final year=10 of  burn-in phase as global environment 
# for each simulation replicate. This will be loaded for other scenarios

cat(paste0("saving burn-in for REP ", REP, "\n"))
save.image(paste0("./data/burnin_",REP,".rda"))
