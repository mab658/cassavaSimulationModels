# Define a function splitPop
# to split genetic resources for two mega-environments

splitPop <- function(pop){
  sample_size = floor(0.5*nInd(pop))
  set.seed(777)
  
  # randomly split data in r
  picked = sample(seq_len(nInd(pop)),size = sample_size)
  ME1 <- pop[picked,] # individuals assigned to Mega Env1
  ME2 <- pop[-picked,]# individuals assigned to Mega Env2
  return(list(ME1=ME1,ME2=ME2))
}