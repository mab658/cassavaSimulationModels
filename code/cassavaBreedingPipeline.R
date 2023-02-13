# Fill/populate  cassava breeding program pipeline with unique individuals from founder population
# This is required to initiate a breeding program by populating each
# breeding stage with genotypes from the same generation.

# In a real breeding program, there are generation difference
# This generation difference will be created later in future testing phase

# total genetic variance remains unchanged because of using the same
# founder parents in each cycle year until burn-in phase where the
# parents are recycle

# year 1  crossing
# year 2 - seedling
# year 3 - clonal evaluation
# year 4 - Preliminary Yield Trial
# year 5 - Advance Yield Trial
# year 6 - Uniform Yield Trial
# year 7 - Variety release

cat("Fill cassava breeding scheme pipeline", "\n")


for (year in 1:7){

  cat("CassavaBreedingPipeline year:",year,"of 7\n")

  # Year 1 - generate F1 in the crossing block
  # make nCrosses=200 bi-parental crosses with
  # nProgeny=25 progeny each to generate F1 population

  F1 <- randCross(pop=parents, nCrosses=nCrosses,
                  nProgeny=nProgeny, simParam=SP)

  # Year 1: Beginning of breeding cycle-
  # Seedling nursery evaluates  the seeds (F1) population
 
  if (year < 7){
    SDN <- setPheno(pop=F1, varE=errVarSDN,reps=repSDN,simParam=SP)
  }

  # year 2  Stage 1 - clonal evaluation trial (CET)
  if (year < 6){
    CET <- selectInd(pop=SDN, nInd=nCET, use="pheno",simParam=SP)
    CET <-  setPheno(pop=CET, varE=errVarCET, reps=repCET,simParam=SP)
  }

  # year 3 - Preliminary Yield Trial (PYT)
  if (year < 5){
    PYT <- selectInd(pop=CET,nInd=nPYT,use="pheno", simParam=SP)
    PYT <-  setPheno(pop=PYT,varE=errVarPYT,reps=repPYT, simParam=SP)
  }

  # year 4 - Advanced Yield Trial (AYT)
  if (year < 4){
    AYT <-  selectInd(pop=PYT, nInd=nAYT, use="pheno", simParam=SP)
    AYT <-   setPheno(pop=AYT,varE=errVarAYT, reps=repAYT,simParam=SP)
  }

  #  year 5 - Uniform Yield Trial I (UYT1)
  if (year < 3){
    UYT <- selectInd(pop=AYT,nInd=nUYT, use="pheno", simParam=SP)
    UYT <- setPheno(pop=UYT, varE=errVarUYT, reps=repUYT,simParam=SP)
  }

  #  - variety release
  if (year < 2){
    # selecting variety for release
    variety <- selectInd(pop=UYT,nInd=nVarietySel, use="pheno",simParam=SP)
  }
} # end loop breeding cycle of cassava