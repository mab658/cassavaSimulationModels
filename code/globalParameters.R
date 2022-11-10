# Define global parameters required for simulation running

nSimRun <- 50 #Number of replications of simulation
burninYears <- 10 # length of common burn-in phase
futureYears <- 10 # length of testing period for future phase
nCycles <- burninYears+futureYears # Number of breeding cycles

# Define parameter for the crossing block
nCrosses <- 200 # number of crosses to be 200
nProgeny <- 50 # number of progeny per cross
famSize <- 5 # selected individuals per family or cross
nVarietySel <- 4 # Number of variety released

# Number of causal sites (QTL) and SNP per each of 18  chrom
nQTL <- 300 # 18 Chrs x 300 QTLs = 5400 sites with effect
nSNP <- 700

# mean and variance Degree of dominance
# contributes to inbreeding depression and dominance variance
ddMean <- 0.20 

# - contribute to dominance variance but not inbreeding depression
ddVar <- 0.10 # Variance dominance degree 


# number of entries per each trial evaluation stages
#nSDN <- 10000
nCET <- 1000
nPYT <- 200
nAYT <- 60
nUYT <- 30


# In AlphaSimR context,number of replicates (rep) is synonymous to the
# number of locations for each evaluation stage
repSDN <- 1
repCET <- 1
repPYT <- 2
repAYT <- 3
repUYT <- 3

### Genomic selection parameters
#Start  training population historical records  of 5 years

startTrainPop <- 8

#Sliding window year limitation => no. of historical records to keep
# after limit year of 3 then remove old records
limityear <- burninYears - (startTrainPop - 1) #do not change

# Estimate mean error from historical data for each yield trial stages
errVarSDN <- 200.5
errVarCET <- 67.6
errVarPYT <- 62.4
errVarAYT <- 41.1
errVarUYT <- 39.5
