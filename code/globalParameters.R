# Define global parameters required for simulation running

burninYears <- 20 # length of common burn-in phase
futureYears <- 10 # length of testing period for future phase
nCycles <- burninYears+futureYears # Number of breeding cycles

# Define parameter for the crossing block
nParents <- 25
nCrosses <- 200 # number of crosses to be 200
nProgeny <- 50 # number of progeny per cross
famSize <- 5 # selected 5individuals per family or cross
nVarietySel <- 4 # Number of variety released


# Number of causal sites (QTL) and SNP per each of 18  chrom
nQTL <- 150 # 18 Chrs x 200 QTLs = 3,600  QTLs causal variants
nSNP <- 500  # simulate 600 x 18 = 10,800  SNPs chip with 10,800  markers


# set the mean and variance of genetic values for additive effect
genMean <- 0.0
varGen <- 1.0

# Set mean (contribute to inbreeding depression)  and variance Degree of dominance for dominance effect
# contributes to inbreeding depression and dominance variance

ddMean <- 0.20 # mean dominance degree 
ddVar <- 0.10 # Variance dominance degree

# set variance of GxE for trait GxE effect
#varGxE <- 0
varGxE <- 10 #(6.5/7.6) * (1/0.438)


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
# The error variance  at SDN is estimated under the assumption that H2 is 5% result to 19
# based on additive genetic variance of 1 for the founder population
# Error from CET from empirical data = 67.6 follow assumption that its H2 is in range of 0.20
# with additive variance of the founder 1, then the error variance is 4.0
# The errors at PYT, AYT and UYT were scaled by 4/67.6

errVarSDN <- 19.0 # H2 = 0.05
errVarCET <- 4.0  # H2 = 0.20 # 67.6
errVarPYT <- 3.7  # (4/67.6* 62.4)
errVarAYT <- 2.4 # (4/67.6 * 41.1)
errVarUYT <- 2.3 # (4/67.6 * 39.5)
