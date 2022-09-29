# # Alternative breeding scenario 2

# This script implemented narrow adaptation breeding program
# by breeding for two separate breeding programs

megaEnv <- c("ME1","ME2") # vector of two mega-environments

splitTrainPop <- splitGen(trainPop) # split training population into 2 for MEs

#SDN <- F1
splitSDN <- splitGen(SDN)# split the SDN population size 10,000  into 2

# split the CET population size 1000 for 2 mega-environments
splitCET <- splitGen(CET) # split the CET population size 1000 for 2 MEs

PYTme1 <- PYT
PYTme2 <- PYT

AYTme1 <- AYT
AYTme2 <- AYT

# loop through each of two mega-environments
for (M in megaEnv){
  if (M=="ME1"){
    SDN <- splitSDN$me1
    trainPop <- splitTrainPop$me1
    gsModel <- RRBLUP(pop=trainPop, traits=1, simParam=SP)

    # Estimate EBV for CET, PYT, and AYT for ME1
    CET <- setEBV(pop=splitCET$me1, solution=gsModel)
    PYT <- setEBV(pop=PYTme1, solution=gsModel)
    AYT <-  setEBV(pop=AYTme1, solution=gsModel)
    source("../code/megaEnv1.R") # simulate for mega-environment 1

  }else { # beginning of mega-environment 2
    SDN <- splitSDN$me2
    trainPop <- splitTrainPop$me2
    gsModel <- RRBLUP(pop=trainPop, traits=1, simParam=SP)

    # Estimate EBV for CET, PYT, and AYT for ME2
    CET <- setEBV(pop=splitCET$me2, solution=gsModel)
    PYT <- setEBV(pop=PYTme2, solution=gsModel)
    AYT <-  setEBV(pop=AYTme2, solution=gsModel)
    source("../code/megaEnv2.R")
  } # end else
}# end for
