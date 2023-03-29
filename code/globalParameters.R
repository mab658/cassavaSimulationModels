# Define global parameters required for simulation running

nSimRun <- 25  # Number of simulation repetitions
burninYears <- 10 # length of common burn-in phase
futureYears <- 10 # length of testing period for future phase
nCycles <- burninYears+futureYears # Number of breeding cycles

# Define parameter for the crossing block
nParents <- 25
nCrosses <- 200 # number of crosses to be 200
nProgeny <- 50 # number of progeny per cross
famSize <- 10 # selected individuals per family or cross
nVarietySel <- 4 # Number of variety released

# Number of causal sites (QTL) and SNP per each of 18  chrom
nQTL <- 150 # 18 Chrs x 150 QTLs = 2700  QTLs causal variants
nSNP <- 600 # simulate 300 x 18 = 5400  SNPs chip with 5400  markers


# set the mean and variance of genetic values for additive effect
genMean <- 0
varGen <- 7.69
#varGen <- 54.7
# Set mean (contribute to inbreeding depression)  and variance Degree of dominance for dominance effect
# contributes to inbreeding depression and dominance variance
ddMean <- 0.20 # mean dominance degree 
ddVar <- 0.10 # Variance dominance degree 


# set variance of GxE for trait GxE effect
#varGxE <- 0
varGxE <- 6.53
#varGxE <- 54.7

#varGxY <- 5 # Genotype by year variance
#varGxL <- 5 #Genotype by location variance
#varGxE <- varGxY + varGxL # Genotype by environment variance


# number of entries per each trial evaluation stages
#nSDN <- 10000
nCET <- 1000
nPYT <- 200
nAYT <- 60
nUYT <- 30

# In AlphaSimR context,number of replicates (rep) is synonymous to the
# number of locations for each evaluation stage

# The effective replication of yield trials to simulate phenotype values
repSDN <- 1
repCET <- 1
repPYT <- 2
repAYT <- 3
repUYT <- 3

# Estimate mean error from historical data for each yield trial stages
# These are the plot-level environmental variance required to simulate phenotypic values
errVarSDN <- 200.5
errVarCET <- 67.6
errVarPYT <- 62.4
errVarAYT <- 41.1
errVarUYT <- 39.5


### Genomic selection parameters
#Start  training population historical records  of 5 years

startTrainPop <- 8

#Sliding window year limitation => no. of historical records to keep
# after limit year of 5 then remove old records
limityear <- burninYears - (startTrainPop - 1) #do not change
