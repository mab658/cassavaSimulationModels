rm(list=ls())
library(AlphaSimR)
library(ASRgenomics)
library(AGHmatrix)
library(asreml)
library(parallel)
library(dplyr)
library(tibble)
library(ggplot2)
library(tidyr)

# setwd("/Users/mab658/Documents/cassavaSimulationModels")
# Note define a cluster using makeCluster for FORK parallel backend do not work 

# set the number of cores and register  parallel backend
# mc.cores <- parallel::detectCores()/2 # number of cores to use

#doParallel::registerDoParallel(cores = nCores) # register parallel backend
mc.cores <- 30

# Load global parameters
source("./code/globalParameters.R")

# load a custom function to split population
# to split genetic resources for two mega-environments
source("./code/splitPop.R")

# Run scripts to load functions  to simulated gxe data  and fit GS model
source("./code/gxeSim.R")
source("./code/gsModelFit.R") # invoke function to fit GBLUP model
source("./code/convProgram.R") # simulate conventional breeding program

# scripts to invoke functions to simulate GS based breeding program
source("./code/broadAdaptation.R")
source("./code/narrowAdaptation.R")


(start.time <- Sys.time())

# Simulate 10-year of burn-in phase for each simulation replicate
for(REP in 1:nSimRun){
  cat("\n Simulation run:",REP,"of", nSimRun,"replication", "\n")

  # Create founder genomes and initial parents population
  source("./code/createParents.R")

  # Fill breeding pipeline with unique individual. from initial parents
  source("./code/cassavaBreedingPipeline.R")

  # Run 10 years of burn-in by advancing years of breeding
  source("./code/burnin.R") # advance evaluation stages by years
}


set.seed(1234, kind = "L'Ecuyer-CMRG")
# Evaluation phase of different breeding scenarios
# mclapply splits these iterations into multiple processes
#  Scenario 1: conventional breeding program

result.conv <- do.call(rbind,mclapply(1:nSimRun, FUN=baseMod, mc.preschedule = TRUE,
                        mc.set.seed = FALSE, mc.cores=mc.cores,mc.silent=FALSE))

distPar.conv <- result.conv %>%
        dplyr::filter(year > burninYears) %>%
        dplyr::select(simRun,year,scenario,nParCET,nParPYT,nParAYT,nParUYT)

 cat("\n Number of trial stages per parental generation \n")
print(distPar.conv)

result.conv <- subset(result.conv,select = -c(nParCET,nParPYT,nParAYT,nParUYT))
print(result.conv)

#write.csv(result.conv, file="conv.csv",row.names=FALSE)


# Scenario 2: Broad adaptation breeding program

set.seed(1234, kind = "L'Ecuyer-CMRG")
result.broad <- do.call(rbind,mclapply(1:nSimRun, FUN=gsBroad, 
			mc.preschedule = TRUE,mc.set.seed = TRUE, 
			mc.cores=mc.cores, mc.silent=FALSE))


# number of trial stages per parental generation
cat("\n Number of trial stages per parental generation \n")

distPar.BA <- result.broad %>%
        dplyr::filter(year > burninYears) %>%
        dplyr::select(simRun,year,scenario,nParCET,nParPYT,nParAYT,nParUYT)

 cat("\n Number of trial stages per parental generation \n")
print(distPar.BA)

result.broad <- subset(result.broad,select = -c(nParCET,nParPYT,nParAYT,nParUYT))
print(result.broad)

#write.csv(result.broad, file="broad.csv",row.names=FALSE)


# Scenario 3: Narrow adaptation - megaEnv1

# This script implemented narrow adaptation breeding program
# with GS to estimate BV in CET based on genomic data to boost
# the selection accuracy at this early stage where phenotype
# information is limited
# GS is used to advance from CET, PYT, AYT, and UYT

set.seed(1234, kind = "L'Ecuyer-CMRG")
result.me1 <- do.call(rbind,mclapply(1:nSimRun,
                    FUN=gsNarrow,ME="ME1",pval1UYT=0.6,pval2UYT=0.9,
                    pval1AYT=0.6,pval2AYT=0.7,
                    pvalPYT=0.7,
                    mc.preschedule = TRUE,mc.set.seed = TRUE, 
                    mc.cores=mc.cores, mc.silent=FALSE))

distPar.me1 <- result.me1 %>%
        dplyr::filter(year > burninYears) %>%
        dplyr::select(simRun,year,scenario,nParCET,nParPYT,nParAYT,nParUYT)

 cat("\n Number of trial stages per parental generation \n")
print(distPar.me1)


cat("\n Simulated data for mega environment 1 - ME1 \n")

result.me1 <- subset(result.me1,select = -c(nParCET,nParPYT,nParAYT,nParUYT))
print(result.me1)



# Scenario 3: Narrow adaptation - megaEnv2
# This script implemented narrow adaptation breeding program for ME2
# with GS to estimate BV in CET based on genomic data to boost
# the selection accuracy at this early stage where phenotype
# information is limited
# GS is used to advance individuals from CET, PYT, AYT, and UYT


cat("\n Modeling scenario 3 narrow adaptation - ME2 enabled GS \n")

set.seed(1234, kind = "L'Ecuyer-CMRG")
result.me2 <- do.call(rbind,mclapply(1:nSimRun,
                    FUN=gsNarrow,ME="ME2",pval1UYT=0.1,pval2UYT=0.4,
                    pval1AYT=0.3,pval2AYT=0.4,
                    pvalPYT=0.3,
                    mc.preschedule=TRUE,mc.set.seed=TRUE, 
                    mc.cores=mc.cores,mc.silent=FALSE))


distPar.me2 <- result.me2 %>%
        dplyr::filter(year > burninYears) %>%
        dplyr::select(simRun,year,scenario,nParCET,nParPYT,nParAYT,nParUYT)

cat("\n # of individuals from evaluation stages that constitute parental candidate - ME2 \n")
#print(distPar.me2)

result.me2 <- subset(result.me2,select = -c(nParCET,nParPYT,nParAYT,nParUYT))


cat("\n Simulated data for mega environment 2 - ME2 \n")
print(result.me2)


######## Saving the simulated data ############

simData <- rbind(result.conv,result.broad,result.me1,result.me2)
write.csv(simData,file=paste0("./data/simData",".csv"),row.names=FALSE)
saveRDS(simData,file=paste0("./data/simData",".rds"))

# rescale cycle year so that last burnin year is zero
simData$year <- -(burninYears-1):burninYears

# compute genetic gain and genetic variance by cycle from last year of burn-in to last future year

simSumm <- simData %>%
  dplyr::filter(year>=0) %>%
  group_by(scenario,year) %>%
  dplyr::summarise(
    genGain = mean(meanGV),
    seGain = sd(meanGV)/sqrt(length(simRun)),
    genVar = mean(varGV),
    seGenVar = sd(varGV)/sqrt(length(simRun)),
    selAccur= mean(accPYT),
    seAccur = sd(accPYT)/sqrt(length(simRun))
)

simSumm$scenario <- factor(simSumm$scenario,
				levels = c("Conv","Broad adaptation","Narrow adaptation"),
				labels = c("Conventional","Broad adaptation","Narrow adaptation"))

cat("\n Summary of genetic gain & variance from last year of burnin to till last cycle \n")
print(simSumm)

# plot of genetic gain

genGainPlot <- ggplot(data=simSumm,aes(x=year,y=genGain,color=scenario))+
  geom_point()+
  geom_line(size=1)+
  geom_errorbar(aes(ymin=genGain-seGain, ymax=genGain+seGain), width=0.2,
                position=position_dodge(0.05))+
  guides(scale="none")+
  theme_bw()+
  theme(legend.position = c(0.03,0.96),
        legend.justification = c("left","top"))+
  scale_x_continuous("Year",limits=c(0,futureYears))+
  scale_y_continuous("Genetic gain",limits=c(0,NA))

# print(genGainPlot)

# save the plot to a file
ggsave("./output/geneticGain.jpeg",height=4.2, width=6.5, units="in", dpi=300)


# plot of genetic variance

genVarPlot <- ggplot(data=simSumm,aes(x=year,y=genVar,color=scenario))+ geom_point()+
  geom_line(size=1)+
  geom_errorbar(aes(ymin=genVar-seGenVar, ymax=genVar+seGenVar), width=0.2,
                position=position_dodge(0.05))+ 
  guides(scale="none")+
  theme_bw()+
  theme(legend.position = c(0.03,0.05),
        legend.justification = c("left","bottom"))+
  scale_x_continuous("Year",limits=c(0,futureYears))+
  scale_y_continuous("Genetic variance",limits=c(0,NA))

# print(genVarPlot)

# save the plot to a file
ggsave("./output/geneticVariance.jpeg",height=4.2, width=6.5, units="in", dpi=300)



# plot distribution of the number of individuals that constitute the parental candidate
parDist <- rbind(distPar.conv,distPar.BA, distPar.me1, distPar.me2)

write.csv(parDist,file=paste0("./data/parentDist",".csv"),row.names=FALSE)
saveRDS(parDist,file=paste0("./data/parDist",".rds"))


parCount <- parDist %>%
  gather(key = stage,value = nPar,-c(simRun,year,scenario))

parCount$scenario <- factor(parCount$scenario,
                           levels = c("Conv","Broad adaptation",
                                      "Narrow adaptation"),
                        labels=c("Conventional","Broad adaptation",
                                 "Narrow adaptation"))

parCount$stage <- factor(parCount$stage,
                        levels = c("nParCET","nParPYT","nParAYT","nParUYT"),
                        labels=c("CET","PYT","AYT","UYT"))

write.csv(parCount,file=paste0("./data/distParentCount",".csv"),row.names=FALSE)

print(parCount)

distParCand  <- ggplot(data = parCount, aes(x=scenario,y=nPar, fill=stage)) +
    geom_boxplot( stat = "boxplot",  outlier.colour = "red", outlier.shape=16,
        outlier.size = 0.5, na.rm=TRUE) +
  labs(x= "Breeding strategy", y= "Parental lines") + theme_bw() +
  theme(axis.title = element_text(colour="black", size=12),
         plot.title = element_text(hjust = 0.5,lineheight=.5,colour="black", size=12),
         axis.text = element_text(face="bold", size=7))

# print(distParCand)
# save the plot to a file
ggsave("./output/distParentalCandidate.jpeg",height=4.2, width=6.5, units="in", dpi=300)



# boxplot of comparing the distribution of  accuracy of selection criterion acros the stages

selAccur <- simData %>%
  dplyr::filter(year>0) %>%
dplyr::select(simRun, year,scenario,accCET,accPYT,accAYT,accUYT) %>%
  gather(key = stage,value = accur,-c(simRun,year,scenario))

selAccur$scenario <- factor(selAccur$scenario,
                           levels = c("Conv","Broad adaptation",
                                      "Narrow adaptation"),
                        labels=c("Conventional","Broad adaptation",
                                 "Narrow adaptation"))

selAccur$stage <- factor(selAccur$stage,
                        levels = c("accCET","accPYT","accAYT","accUYT"),
                        labels=c("CET","PYT","AYT","UYT"))

distSelAccur <- ggplot(data = selAccur, aes(x=scenario,y=accur, fill=stage)) +
    geom_boxplot( stat = "boxplot",  outlier.colour = "red", outlier.shape=16, 
	outlier.size = 0.5, na.rm=TRUE) +
  labs(x= "Breeding strategy", y= "Selection accuracy") + theme_bw() +
  theme(axis.title = element_text(colour="black", size=12),
         plot.title = element_text(hjust = 0.5,lineheight=.5,colour="black", size=12),
         axis.text = element_text(face="bold", size=7))
   
# print(distSelAccur)
# save the plot to a file
ggsave("./output/distSelectionAccuracy.jpeg",height=4.2, width=6.5, units="in", dpi=300)


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

print(selAccur)

selAccPlot <- ggplot(data=selAccur,aes(x=scenario,y=accur,fill=stage)) + 
geom_bar(stat = "identity", position = "dodge")+
  geom_text(aes(x=scenario,y=accur,label=round(accur,2)),
            position = position_dodge(width = 0.9),vjust=-0.25)+
  labs(x="Breeding programme", y="Selection accuracy")
 
#print(selAccPlot)
 
 # save the plot to a file
ggsave("./output/selAccuracyMean.jpeg",height=4.2, width=6.5, units="in", dpi=300)




# plot line trends of the accuracy at PYT stage

selAccurPlot <- ggplot(data=simSumm,aes(x=year,y=selAccur,color=scenario))+
  geom_point()+
  geom_line(size=1)+
  geom_errorbar(aes(ymin=selAccur-seAccur, ymax=selAccur+seAccur), width=0.2,
                position=position_dodge(0.05))+
  guides(scale="none")+
  theme_bw()+
  theme(legend.position = c(0.03,0.96),
        legend.justification = c("left","top"))+
  scale_x_continuous("Year",limits=c(0,futureYears))+
  scale_y_continuous("Selection accuracy",limits=c(0,NA))

# print(selAccurPlot)

# save the plot to a file
ggsave("./output/predAccuracy.jpeg",height=4.2, width=6.5, units="in", dpi=300)




# fit linear regression of genetic gain on the breeding cycle of evaluation phase
# to calculate rate of genetic gain per year

fitConv <- lm(genGain ~ year, data = simSumm[simSumm$scenario=='Conventional' & simSumm$year > 0,])
cat("Rate of genetic gain per year for conventional breeding program is",
    round(coef(fitConv)[2],2),"\n")


fitBroad <- lm(genGain ~ year, data = simSumm[simSumm$scenario=='Broad adaptation' & simSumm$year > 0,])
cat("Rate of genetic gain per year for broad adaptation breeding program is",
    round(coef(fitBroad)[2],2),"\n")


fitNarrow <- lm(genGain ~ year, data = simSumm[simSumm$scenario=='Narrow adaptation' & simSumm$year > 0,])
cat("Rate of genetic gain per year for narrow adaptation breeding program is",
    round(coef(fitNarrow)[2],2),"\n")

(end.time <- Sys.time())
(end.time - start.time)

doParallel::stopImplicitCluster()
