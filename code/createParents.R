# Simulate founder genomes and define trait genetic architecture
# from coalescent simulator called runMac2
# i.e. Create haplotypes sequence genomes for founder population of 50 outbred
# individuals to mimic crop evolution

founderGenomes <- runMacs2(nInd=nParents,
                       nChr=18,
                       segSites=nQTL+nSNP,
                       Ne = 100,
                       #bp=1e+08, # physical length (bp)
                       #genLen=1.43, # genetic length (M)
                       #mutRate = 2.5e-08,
                       inbred= FALSE,
                       #species = "GENERIC",
                       ploidy = 2L
)

# genLen = genetic length of chromosomes in Morgan e.g 1.43 M
# bp = SNPs Physical position or length on the chromosomes of 1 x 108 base pairs

# Create an object (SP) to set the global simulation parameters
# for founder haplotypes sequences genomes
# this set the genotype to phenotype mapping

SP <- SimParam$new(founderGenomes)


# define a trait with additive, dominance, and GxE genetic architecture
# mean = mean of additive genetic value
# var = variance of genetic value

# Degree of dominance: 0 = no dominance, between (0,1) = partial dominance,
# 1= complete dominance, > 1 = overdominance
# The average degree of dominance is estimated as the square root of the average
# squared degree of dominance.

#meanDD - mean dominance degree > 0 means trait with directional dominance 0 all genes
# show dominance in the same direction
#varDD - Variance of dominance degree

# we also modelled both environmental variance (varEnv=0) and
# GxE interaction variance (varGxE=15 or 30)

# Sets restrictions on which segregating sites
# can serve as SNP or QTL
# `restrSegSites prevents SNP from also being QTL

SP$restrSegSites(nQTL, nSNP, overlap=FALSE) # maxQTL=nQTL and maxSNp=nSNP


SP$addTraitADG(nQtlPerChr=nQTL, 
		mean=genMean, var = varGen,
              	meanDD=ddMean,varDD=ddVar,
              	varGxE=varGxE)


# set a phenotype for the initial parents in the founder population
#SP$setVarE(H2=0.30)

#  Phenotype the initial parental population created from founder genomes
#  parents = setPheno(parents,  varE = VarE)

# add SNPs chip for all simulated SNPS markers in the population
# SNPs are generally described by the number of SNP position the assay
#creates a marker set with nSNP  per chromosome

SP$addSnpChip(nSnpPerChr=nSNP)

#SP$setTrackRec(TRUE) # keep track records


# create initial parental base population from founder genomes
# from which unique individuals in each evaluation stage are generated
# All the parents are effective the same generation being from same founder genomes

# beginning of modeling the breeding program.Create  parental population from founder genomes
parents <- newPop(founderGenomes, simParam=SP)
