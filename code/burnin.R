# Load global parameters
source("./code/globalParameters.R")

# Create founder genomes and initial parents population
source("./code/createParents.R")

# Fill breeding pipeline with unique individual from initial parents
source("./code/cassavaBreedingPipeline.R")

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
                     "phenoNA","genoNA")


# we track the progress of breeding program
# by initializing variables to save the mean and variance of genetic values
simRun <- vector("numeric",length=nCycles)
meanGV <- vector("numeric",length=nCycles)
varGV <- vector("numeric",length=nCycles)

h2 <- H2 <- vector("numeric",length=nCycles)
ebv <- egv <- vector("numeric",length=nCycles)

# placeholder for selection accuracy at each stage
accSDN <- accCET <- accPYT <- accAYT <- accUYT <- vector("numeric",length=nCycles)


# pipeline  to avoid copying data

# track number of individuals selected as parent in each trial stage per year
nParCET <- nParPYT <- nParAYT <- nParUYT <- vector("numeric",length=nCycles)


# create a placeholder to store p-value range to denote environmental covariate of GxE  for each location
allLocPvals <- numeric(9)
#pvalRange <- as.data.frame(matrix(c(0.1, 0.5, 0.15, 0.55, 0.2, 0.6, 0.2, 0.6,
#                          0.3, 0.7, 0.4, 0.8, 0.4, 0.8, 0.45, 0.85, 0.5, 0.9),
#                          ncol=2,byrow = TRUE))


eCovRange <- as.data.frame(matrix(c(-1.9, -0.7,-0.9,0.3,-0.9,0.3,-0.8,0.4,-0.6,0.6,-0.3,0.9,
				0.0,1.2,0.1,1.3,0.9,2.1),
				ncol=2,byrow = TRUE))



# Run 10 years of burn-in (Cycle years) as a common starting point
# for different scenarios.
# This marks the beginning of a new breeding cycle where new parents are
# formed and taking to crossing block

# Advance breeding program by Working backwards through pipeline
# to avoid copying data

for(year in 1:burninYears){#Change to any number of desired years
	cat("Burnin phase cycle year:",year,"of",burninYears, "\n")


	# recycle new parents in each year of burn-in phase
        # by selecting best individuals from joint (UYT and AYT) population

        parents <- selectInd(pop = c(UYT,AYT), nInd=nParents, trait = 1,use="pheno",simParam=SP)

  	# randomly sample from environmental covariate and convert to p-value 
	# to be assigned to each location in a cycle year
  	for (loc in 1:9){
       		allLocPvals[loc] <- pnorm(runif(1, eCovRange[loc,1], eCovRange[loc,2]))
  	}

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
	CETrec <- gxeSim(pvalVec=allLocPvals[5], pop = CET,
                        varE = errVarCET,nreps = repCET, locID = 5,year=year)
        CET <- CETrec[[1]]


  	# Seedling Nursery (SDN)
	SDN <- F1
	SDNrec  <- gxeSim(pvalVec=allLocPvals[5], pop = SDN,
                         varE = errVarSDN,nreps = repSDN, locID = 5,year=year)
        SDN <- SDNrec[[1]]

  	# recycle new parents in each year of burn-in phase
  	# by selecting best individuals from joint (UYT and AYT) population

  	# crossing block
  	F1 <- randCross(pop = parents, nCrosses = nCrosses,
                  nProgeny = nProgeny, simParam = SP)

  	# report the mean and variance of genetic values at PYT in each year
  
	meanGV[year] <- meanG(PYT)
  	varGV[year] <- varG(PYT)
 
} # end 10-year burn-in

# Save the state of simulation at final year=10 of  burn-in phase as global environment 
# for each simulation replicate. This will be loaded for other scenarios

cat(paste0("saving burn-in for REP ", REP, "\n"))
save.image(paste0("./data/burnin_",REP,".rda"))

