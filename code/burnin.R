# After populating the cassava breeding pipeline
# we implement the burn-in phase by advancing population of individuals
# through the breeding scheme over years where the parents are recycled.
# Burnin phase serves as common starting point for
# different breeding scenarios (programs)
# Advance breeding by year by working backwards through
# pipeline  to avoid copying data

# create a placeholder of training records for genomic selection (GS)
#trainPop <- data.frame()
trainRec <- vector(mode="list",length=8)
names(trainRec) <- c("phenoBA","genoBA",
                     "phenoME1","genoME1",
                     "phenoME2","genoME2",
                     "phenoNA","genoNA")


# we track the progress of breeding program
# by initializing variables to save the mean and variance of genetic values
simRun <- vector("numeric",length=nCycles)
meanGV <- vector("numeric",length=nCycles)
varGV <- vector("numeric",length=nCycles)
eGain <- vector("numeric",length=nCycles)

H2 <- vector("numeric",length=nCycles)

# placeholder for selection accuracy at each stage
accCET <- accPYT <- accAYT <- accUYT <- vector("numeric",length=nCycles)


# pipeline  to avoid copying data

# track number of individuals selected as parent in each trial stage per year
nParCET <- nParPYT <- nParAYT <- nParUYT <- vector("numeric",length=nCycles)




# Run 10 years of burn-in (Cycle years) as a common starting point
# for different scenarios.
# This marks the beginning of a new breeding cycle where new parents are
# formed and taking to crossing block

# Advance breeding program by Working backwards through pipeline
# to avoid copying data

# Generate or replace presample p-value for each cycle year 
# in the burn-in and future evaluation steps
 
allLocPvals <- numeric(9)
pvalRange <- as.data.frame(matrix(c(0.1, 0.5, 0.15, 0.55, 0.2, 0.62, 0.2, 0.68,
                          0.3, 0.7, 0.4, 0.8, 0.4, 0.8, 0.45, 0.85, 0.5, 0.9),
                          ncol=2,byrow = TRUE))

# Assign p-value to each location in a cycle year
for (loc in 1:9){
       allLocPvals[loc] <- runif(1, pvalRange[loc, 1], pvalRange[loc, 2])
 }



for(year in 1:burninYears){#Change to any number of desired years
  cat("Burnin phase cycle year:",year,"of",burninYears, "\n")
  
  # Select variety to be released
  variety <- selectInd(pop = UYT, nInd = nVarietySel,trait = 1, use = "pheno",simParam = SP)

  # Uniform Yield Trial (UYT)

  UYT <- selectInd(pop = AYT, nInd = nUYT, trait = 1, use = "pheno",simParam = SP)
 
  UYTrec <- gxeSim(pvalVec=allLocPvals[-5], pop = UYT,
                varE = errVarUYT, nreps = repUYT,locID=c(1:4,6:9),year=year)
  UYT <- UYTrec[[1]]


  # Advance Yield Trial (AYT)
  AYT <- selectInd(pop = PYT, nInd = nAYT, trait = 1, use = "pheno",simParam = SP)
 

  # Invoke the function to phenotype selected AYT clones in 4 locations
  AYTrec <- gxeSim(pvalVec=allLocPvals[c(3,4,6,7)],pop = AYT,
                varE = errVarAYT, nreps = repAYT, locID = c(3,4,6,7),year=year)
  AYT <- AYTrec[[1]]



  # Preliminary  Yield  Trial (PYT)
  PYT <- selectInd(pop = CET, nInd = nPYT, trait = 1, use = "pheno", simParam = SP)

   # Invoke the function to phenotype selected PYT clones in 2 locations
  PYTrec <- gxeSim(pvalVec=allLocPvals[c(3,7)],pop = PYT,
                varE = errVarPYT, nreps = repPYT, locID=c(3,7),year=year)
  PYT <- PYTrec[[1]]


  #  Clonal Evaluation Trial (CET)
  CET <- selectWithinFam(pop = SDN, nInd = famSize, trait = 1, use = "pheno", simParam = SP)
  CET <- selectInd(pop = CET, nInd = nCET,trait = 1, use = "pheno", simParam = SP)
  #Invoke the function to phenotype selected individuals in 1 locations
  CETrec <- gxeSim(pvalVec=allLocPvals[5], pop = CET,
                   varE = errVarCET,nreps = repCET,locID=5,year=year)

  CET <- CETrec[[1]]


  # Seedling Nursery (SDN)
  SDNrec  <- gxeSim(pvalVec=allLocPvals[5], pop = F1,
                   varE = errVarSDN,nreps = repSDN, locID = 5,year=year)
  
  SDN <- SDNrec[[1]]



  # recycle new parents in each year of burn-in phase
  # by selecting best individuals from joint (UYT and AYT) population

  parents <- selectInd(pop = c(UYT,AYT), nInd=nParents, trait = 1,
                       use="pheno",simParam=SP)

  # crossing block
  F1 <- randCross(pop = parents, nCrosses = nCrosses,
                  nProgeny = nProgeny, simParam = SP)

  
  # report the mean and variance of genetic values at PYT in each year
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
   # splitCET <- splitPop(CET)
    
    # collect genotype & phenotype of inds across stages for ME1
    #CETrecME1 <-  CETrec[[2]][CETrec[[2]]$gen %in% splitCET$ME1@id,]

    trainRec$phenoME1 <- rbind(CETrec[[2]],PYTrec[[2]][PYTrec[[2]]$loc=="L7",],
                               AYTrec[[2]][AYTrec[[2]]$loc %in% c("L6","L7"),],
                               UYTrec[[2]][UYTrec[[2]]$loc %in% c("L6","L7","L8","L9"),])

     

    # genotyping data at beginning of TP for ME1    

    trainRec$genoME1 <- pullSnpGeno(pop = c(CET, PYT, AYT, UYT), simParam = SP)
    
       
    # collect phenotype and genotype data of individuals for ME2 
    
    #CETrecME2 <-  CETrec[[2]][CETrec[[2]]$gen %in% splitCET$ME2@id,]
    #trainRec$phenoME2 <- rbind(CETrecME2,PYTrec[[2]],
    #                           AYTrec[[2]],UYTrec[[2]])

    trainRec$phenoME2 <- rbind(CETrec[[2]],PYTrec[[2]][PYTrec[[2]]$loc=="L3",],
                               AYTrec[[2]][AYTrec[[2]]$loc %in% c("L3","L4"),],
                               UYTrec[[2]][UYTrec[[2]]$loc %in% c("L1","L2","L3","L4"),])
   

    # genotyping data at beginning of TP for ME2
    trainRec$genoME2 <- pullSnpGeno(pop=c(CET, PYT, AYT, UYT),simParam = SP)
    
    
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
    # splitCET <- splitPop(CET) # Split the CET to constitute TP for ME1 and ME2

    #CETrecME1 <-  CETrec[[2]][CETrec[[2]]$gen %in% splitCET$ME1@id,]

    trainRec$phenoME1 <- rbind(trainRec$phenoME1,
                               CETrec[[2]],PYTrec[[2]][PYTrec[[2]]$loc=="L7",],
                               AYTrec[[2]][AYTrec[[2]]$loc %in% c("L6","L7"),],
                               UYTrec[[2]][UYTrec[[2]]$loc %in% c("L6","L7","L8","L9"),])


    trainRec$genoME1 <- rbind(trainRec$genoME1,pullSnpGeno(pop=c(CET, PYT, AYT, UYT),simParam=SP))
    
    # update TP for ME2
    
    #CETrecME2 <-  CETrec[[2]][CETrec[[2]]$gen %in% splitCET$ME2@id,]

    trainRec$phenoME2 <- rbind(trainRec$phenoME2,
                               CETrec[[2]],PYTrec[[2]][PYTrec[[2]]$loc=="L3",],
                               AYTrec[[2]][AYTrec[[2]]$loc %in% c("L3","L4"),],
                               UYTrec[[2]][UYTrec[[2]]$loc %in% c("L1","L2","L3","L4"),])
    
    trainRec$genoME2 <- rbind(trainRec$genoME2,  pullSnpGeno(pop=c(CET, PYT, AYT, UYT),simParam=SP))

    
  } # end of elseif
} # end 10-year burn-in
 
write.csv(trainRec$phenoBA,file="trainRecPhenoBA.csv", row.names=FALSE)
write.csv(trainRec$phenoME1,file="trainRecPhenoME1.csv", row.names=FALSE)
write.csv(trainRec$phenoME2,file="trainRecPhenoME2.csv", row.names=FALSE)

# Save the state of simulation at final year=10 of  burn-in phase as global environment 
# for each simulation replicate. This will be loaded for other scenarios

cat(paste0("saving burn-in for REP ", REP, "\n"))
save.image(paste0("./data/burnin_",REP,".rda"))
