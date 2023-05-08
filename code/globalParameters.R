# Define global parameters required for simulation running

burninYears <- 20 # length of common burn-in phase
futureYears <- 10 # length of testing period for future phase
nCycles <- burninYears+futureYears # Number of breeding cycles

# Define parameter for the crossing block
nParents <- 25
nCrosses <- 200 # number of crosses to be 200
nProgeny <- 50 # number of progeny per cross
famSize <- 5 # selected 20 individuals per family or cross
nVarietySel <- 4 # Number of variety released


# Number of causal sites (QTL) and SNP per each of 18  chrom
nQTL <- 300 # 18 Chrs x 200 QTLs = 3,600  QTLs causal variants
nSNP <- 700 # simulate 600 x 18 = 10,800  SNPs chip with 10,800  markers


# set the mean and variance of genetic values for additive effect
genMean <- 0.0
varGen <- 1

# Set mean (contribute to inbreeding depression)  and variance Degree of dominance for dominance effect
# contributes to inbreeding depression and dominance variance
ddMean <- 0.20 # mean dominance degree 
ddVar <- 0.10 # Variance dominance degree
#ddVar <- 0.50 # Variance dominance degree


# set variance of GxE for trait GxE effect

varGxY <- 1
varGxL <- 1
varGxE <- varGxY + varGxL

#varGxE <- 10
#varGxE <- 46.4 #54.7*6.53/7.69

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
errVarSDN <- 19 # H2=0.05
errVarCET <- 1.08 # H2= 0.2 error=4
errVarPYT <- 1.33 # H2 = 0.43
errVarAYT <- 1.13 # H2 = 0.47
errVarUYT <- 1.38 # H2 = 0.42 


### Genomic selection parameters
#Start  training population historical records  of 3 years

startTrainPop <- 16

#Sliding window year limitation => no. of historical records to keep
limityear <-  burninYears - (startTrainPop - 1) # after limit year, remove old record
