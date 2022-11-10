rm(list=ls())
library(AlphaSimR)
library(rrBLUP)
library(doParallel)
library(dplyr)
library(ggplot2)
library(tidyr)


# set the number of cores to be used

cl <- makeCluster(40, # number of cores to use
                         type = "FORK") # type of cluster
# register doParallel as our parallel backend
registerDoParallel(cores = cl)


# define function to extract gv
gvADG <- function(pop,pval){
	return(setPheno(pop,pval,vE=0))}

# Load global parameters
source("./code/globalParameters.R")

# load a custom function to split population
# to split genetic resources for two mega-environments
source("./code/splitPop.R")

# load a function to generate simulated gxe data
# across multipe locations
source("./code/gxeSim.R")


(start.time <- Sys.time())

# Simulation runs for 10-year burn-in phase for nSimRun replicates
for(REP in 1:nSimRun){  
  cat("\n Simulation run:",REP,"of", nSimRun,"replication", "\n")
  
  # Create founder genomes and initial parents population
  source("./code/createParents.R")
  
  # Fill breeding pipeline with unique individual. from initial parents
  source("./code/cassavaBreedingPipeline.R")
  
  # Run 10 years of burn-in by advancing years of breeding
  source("./code/burnin.R") # advance evaluation stages by years
}


# Evaluation phase of different breeding scenarios


# Scenario 1: conventional breeding program
result.conv <- foreach(REP=icount(nSimRun),
                       .packages=c("AlphaSimR"),
                       .combine=rbind,
                       .multicombine=TRUE,
                       .verbose=TRUE
) %dopar%{

   cat("\n Modeling scenario 1 conventional breeding scheme \n")

   cat("\n Simulation run:",REP,"of", nSimRun,"replication", "\n")

  # load the burnin phase for the second scenario
  load(paste0("./data/burnin_",REP,".rda"))

  # invoke the script to advance year of breeding for genomic selection
  source("./code/conv.R",local=T)
  return(simParms)
}

#print(result.conv)


# Scenario 2: broad adaptation

result.broad <- foreach(REP=icount(nSimRun),
                       .packages=c("AlphaSimR"),
                       .combine=rbind,
                       .multicombine=TRUE,
                       .verbose=TRUE
) %dopar%{

  cat("\n Modeling scenario 2 broad adapatation enabled GS \n")

  cat("\n Simulation run:",REP,"of", nSimRun,"replication", "\n")

  # load the burnin phase for the second scenario
  load(paste0("./data/burnin_",REP,".rda"))

  # invoke the script to advance year of breeding for genomic selection
  source("./code/broadAdaptation.R",local = T)
  return(simParms)
}

#print(result.broad)



# Scenario 3: Narrow adaptation - megaEnv1

# This script implemented narrow  adaptation BP
# with GS to estimate BV in CET based on genomic data to boost 
# the selection accuracy at this early stage where phenotype 
# information is limited
# GS is used to advance from CET, PYT, AYT, and UYT


result.me1 <- foreach(REP=icount(nSimRun),
                       .packages=c("AlphaSimR"),
                       .combine=rbind,
                       .multicombine=TRUE,
                       .verbose=TRUE
) %dopar%{

  cat("\n Modeling scenario 3 narrow  adaptation enabled GS \n")
  
  cat("\n Simulation run:",REP,"of", nSimRun,"replication", "\n")

  # load the burnin phase for the second scenario
  load(paste0("./data/burnin_",REP,".rda"))
  
  # Train genomic selection (GS) model and safe the model fit
  gsModel <- RRBLUP(pop=trainPopME1, traits=1, use="pheno", simParam=SP)
  
  CET = setEBV(pop=CET, solution=gsModel,simParam=SP)
  PYT = setEBV(pop=PYT, solution=gsModel,simParam=SP)
  AYT = setEBV(pop=AYT, solution=gsModel,simParam=SP)
  UYT = setEBV(pop=UYT, solution=gsModel,simParam=SP)
  trainPop <- trainPopME1
  
  source("./code/megaEnv.R",local=T) # run function in local env
  # Invoke megaEnv function to advance narrow-adaptation for megaEnv1
  simParms <- megaEnv(pval1UYT=0.6,pval2UYT=0.9,
                    pval1AYT=0.6,pval2AYT=0.7,
                    pvalPYT=0.7,pvalCET=0.5, 
                    pvalSDN=0.5)
  # write.csv(output,file=paste0("./data/narrowME1_parms","_",REP,".csv"),row.names=FALSE)
  return(simParms)
} # end of doParallel
#print(result.me1)


# Scenario 3: Narrow adaptation - megaEnv2
# This script implemented narrow  adaptation BP
# with GS to estimate BV in CET based on genomic data to boost
# the selection accuracy at this early stage where phenotype
# information is limited
# GS is used to advance from CET, PYT, AYT, and UYT


result.me2 <- foreach(REP=icount(nSimRun),
                       .packages=c("AlphaSimR"),
                       .combine=rbind,
                       .multicombine=TRUE,
                       .verbose=TRUE
) %dopar%{

  cat("\n Modeling scenario 3 narrow  adaptation enabled GS \n")

  cat("\n Simulation run:",REP,"of", nSimRun,"replication", "\n")

  # load the burnin phase for the second scenario
  load(paste0("./data/burnin_",REP,".rda"))

  # Train genomic selection (GS) model and safe the model fit
  gsModel <- RRBLUP(pop=trainPopME2, traits=1, use="pheno", simParam=SP)

  CET = setEBV(pop=CET, solution=gsModel,simParam=SP)
  PYT = setEBV(pop=PYT, solution=gsModel,simParam=SP)
  AYT = setEBV(pop=AYT, solution=gsModel,simParam=SP)
  UYT = setEBV(pop=UYT, solution=gsModel,simParam=SP)
  trainPop <- trainPopME2

  source("./code/megaEnv.R",local=T) # run function in local env

# Invoke megaEnv function to advance narrow-adaptation for megaEnv2
  simParms <- megaEnv(pval1UYT=0.1,pval2UYT=0.4,
                    pval1AYT=0.3,pval2AYT=0.4,
                    pvalPYT=0.3,pvalCET=0.5,
                    pvalSDN=0.5)
 # write.csv(simParms,file=paste0("./data/narrowME2_parms","_",REP,".csv"),row.names=FALSE)
  return(simParms)
}
#print(result.me2)


######## Saving the simulated data ############

simData <- rbind(result.conv,result.broad,result.me1,result.me2)
write.csv(simData,file=paste0("./data/simData",".csv"),row.names=FALSE)
saveRDS(simData,file=paste0("./data/simData",".rds"))

# rescale cycle year so that last burnin year is zero
simData$year <- -(burninYears-1):burninYears

# compute mean genetic value (genetic gain) and genetic variance

simSumm <- simData %>%
  dplyr::filter(year>=0) %>%
  group_by(scenario,year) %>%
  dplyr::summarise(
    gvMean = mean(meanGV),
    sDgvMean = sd(meanGV),
    genGain = mean(gGain),
    genVar = mean(varGV),
    sDgenVar=sd(varGV))

simSumm$scenario <- factor(simSumm$scenario,
				levels = c("Conv","Broad adaptation","Narrow adaptation"),
				labels = c("Conventional","Broad adaptation","Narrow adaptation"))

# plot of genetic gain

gGainPlot <- ggplot(data=simSumm,aes(x=year,y=gvMean,color=scenario))+
  geom_point()+
  geom_line(size=1)+
  geom_errorbar(aes(ymin=gvMean-sDgvMean, ymax=gvMean+sDgvMean), width=.2,
                position=position_dodge(0.05))+
  guides(scale="none")+
  theme_bw()+
  theme(legend.position = c(0.03,0.96),
        legend.justification = c("left","top"))+
  scale_x_continuous("Year",limits=c(0,futureYears))+
  scale_y_continuous("Mean genetic value",limits=c(0,NA))

#print(gGainPlot)

# save the plot to a file
ggsave("./output/geneticGain.jpeg",height=4.2, width=6.5, units="in", dpi=300)


# plot of mean relative genetic variance
varGplot <- ggplot(data=simSumm,aes(x=year,y=genVar,color=scenario))+ geom_point()+
  geom_line(size=1)+
  geom_errorbar(aes(ymin=genVar-sDgenVar, ymax=genVar+sDgenVar), width=.2,
                position=position_dodge(0.05))+ 
  guides(scale="none")+
  theme_bw()+
  theme(legend.position = c(0.03,0.05),
        legend.justification = c("left","bottom"))+
  scale_x_continuous("Year",limits=c(0,futureYears))+
  scale_y_continuous("Mean genetic variance",limits=c(0,NA))
print(varGplot)

# save the plot to a file
ggsave("./output/varGeneticGain.jpeg",height=4.2, width=6.5, units="in", dpi=300)



# boxplot of selection accuracy

selAccur <- simData %>%
  dplyr::filter(year>0) %>%
dplyr::select(simRun, year,scenario,accCET,accPYT,accAYT,accUYT) %>%
  gather(key = stage,value = accur,-c(simRun, year,scenario))

selAccur$scenario <- factor(selAccur$scenario,
                           levels = c("Conv","Broad adaptation",
                                      "Narrow adaptation"),
                        labels=c("Conventional","Broad adaptation",
                                 "Narrow adaptation"))

selAccur$stage <- factor(selAccur$stage,
                        levels = c("accCET","accPYT","accAYT","accUYT"),
                        labels=c("CET","PYT","AYT","UYT"))

ggplot(data = selAccur, aes(x=scenario,y=accur, fill=stage)) +
    geom_boxplot( stat = "boxplot",  outlier.colour = "red", outlier.shape=16, 
	outlier.size = 0.5, na.rm=TRUE) +
  labs(x= "Breeding strategy", y= "Selection accuracy") + theme_bw() +
  theme(axis.title = element_text(colour="black", size=12),
         plot.title = element_text(hjust = 0.5,lineheight=.5,colour="black", size=12),
         axis.text = element_text(face="bold", size=7))
   

# save the plot to a file
ggsave("./output/selectionAccuracy.jpeg",height=4.2, width=6.5, units="in", dpi=300)


# plot Mean of selection accuracy
selAccur <- simData %>%
  dplyr::filter(year>0) %>%
  group_by(scenario) %>%
  dplyr::summarise(
    accCET=mean(accCET),
    accPYT=mean(accPYT),
    accAYT=mean(accAYT),
    accUYT=mean(accUYT)) %>%
  gather(key = stage,value = accur,-c(scenario))

selAccur$scenario <- factor(selAccur$scenario,
                          levels = c("Conv","Broad adaptation","Narrow adaptation"),
                          labels = c("Conventional","Broad adaptation","Narrow adaptation"))

selAccur$stage <- factor(selAccur$stage,
			levels = c("accCET","accPYT","accAYT","accUYT"),
			labels=c("CET","PYT","AYT","UYT"))


selAccPlot <- ggplot(data=selAccur,aes(x=scenario,y=accur,fill=stage)) + 
geom_bar(stat = "identity", position = "dodge")+
  geom_text(aes(x=scenario,y=accur,label=round(accur,2)),
            position = position_dodge(width = 0.9),vjust=-0.25)+
  labs(x="Breeding programme", y="Selection accuracy")
  
 # save the plot to a file
ggsave("./output/selAccuracyMean.jpeg",height=4.2, width=6.5, units="in", dpi=300)


(end.time <- Sys.time())
(end.time - start.time)
stopCluster(cl)
