# Simulate founder genomes and define trait genetic architecture
# from coalescent simulator called runMac2
# i.e. Create haplotypes sequence genomes for founder population of 50 outbred
# individuals to mimic crop evolution

founderGenomes <- runMacs2(nInd=50,
                       nChr=18,
                       segSites=nQTL+nSNP,
                       Ne = 50,
                       bp=8e+08, # physical length (bp)
                       genLen=1.43, # genetic length (M)
                       mutRate = 2e-09,
                       #histNe = c(500, 1500, 6000, 12000, 1e+05),
                       #histGen =  c(100, 1000, 10000, 1e+05, 1e+06),
                       inbred= FALSE,
                       #species = "GENERIC",
                       ploidy = 2L,
                       nThreads = 6
)

# genLen = genetic length of chromosomes in Morgan e.g 1.43 M
# bp = SNPs Physical position or length on the chromosomes of 8 x 108 base pairs

# Create an object (SP) holding the global simulation parameters
# for founder haplotypes sequences genomes
# this set the genotype to phenotype mapping
SP <- SimParam$new(founderGenomes)

# Sets restrictions on which segregating sites
# can serve as SNP or QTL
# `restrSegSites prevents SNP from also being QTL
SP$restrSegSites(nQTL, nSNP, overlap = FALSE) # maxQTL=nQTL and maxSNp=nSNP

# add SNPs chip for all simulated SNPS markers in the population
# SNPs are generally described by the number of SNP position the assay
#creates a marker set with nSNP  per chromosome
# for fitting GS model
SP$addSnpChip(nSnpPerChr=nSNP)

# define a trait with additive, dominance, and GxE genetic architecture
# mean = mean of additive genetic value
# var = variance of genetic value

# Degree of dominance: 0 = no dominance, between (0,1) = partial dominance,
# 1= complete dominance, > 1 = overdominance
#meanDD - mean dominance degree > 0 means trait with directional dominance 0 all genes
# show dominance in the same direction
#varDD - Variance of dominance degree

# we also modelled both environmental variance (varEnv=0) and
# GxE interaction variance (varGxE=0.5)

SP$addTraitADG(nQtlPerChr=nQTL, mean=0, var=1,
              meanDD=ddMean,varDD=ddVar,
              varEnv=0,varGxE=0.5)

SP$setTrackRec(TRUE) # keep track records

# Model the trait with a broad-sense heritability of 0.3
# for evaluation in a single location
#SP$setVarE(H2 = 0.30)

# create initial parental base population from founder genomes
# from which unique individuals in each evaluation stage are generated
# All the parents are effective the same generation
# being from same founder genomes

# begining of modeling the breeding program
parents <- newPop(founderGenomes, simParam=SP)