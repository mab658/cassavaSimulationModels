# Define global parameters required for simulation running

nSimRun <- 5 #Number of replications of simulation
burninYears <- 10 # length of common burn-in phase
futureYears <- 10 # length of testing period for future phase
nCycles <- burninYears+futureYears # No of breeding cycles

# Define parameter for the crossing block
nCrosses <- 200 # number of crosses to be 200
nProgeny <- 50 # number of progeny per cross
famSize <- 5 # selected individuals per family or cross
nVarietySel <- 4 # Number of variety released

# Number of causal sites (QTL) and SNP per each of 18  chrom
nQTL = 300 # 18 chr x 300 QTLs = 5400 alleles with an effect
nSNP = 700

# mean and variance Degree of dominance
ddMean = 0.20 # contributes to inbreeding depression and dominance variance
ddVar = 0.16 # Variance


# number of entries per each trial evaluation stages
#nSDN <- 10000
nCET <- 1000
nPYT <- 100
nAYT <- 40
nUYT <- 30


# In AlphaSimR context,number of replicates (rep) is synoymous to the
# number of locations for each evaluation stage
repSDN <- 1
repCET <- 1
repPYT <- 2
repAYT <- 3
repUYT <- 4


### Genomic selection parameters
#Start  training population historical records  of 3 years
startTrainPop <- 8
#Sliding window year limitation => no. of historical records to keep
# after limit year of say 3 then remove old records
limityear <- burninYears - (startTrainPop - 1) #do not change


# Estimate mean error from historical data for each yield trial stages
errVarSDN <- 200.5
errVarCET <- 67.6
errVarPYT <- 62.4
errVarAYT <- 41.1
errVarUYT <- 39.5
