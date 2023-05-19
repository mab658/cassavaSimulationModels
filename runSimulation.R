# Clean the working environment
rm(list = ls())

# setwd("/Users/mab658/Documents/cassavaSimulationModels")

# list of required packages
packages_used <- c("AlphaSimR", "ASRgenomics", "AGHmatrix", "asreml",
"tidyverse","parallel", "tidyr","tibble","ggplot2","ggpubr")

ip <- installed.packages()
all_packages_installed <- TRUE

for (package in packages_used){
	if (!(package %in% ip[,"Package"])){
    		print(paste("Please install package", package))
    		all_packages_installed <- FALSE
	} else {
		library(package, character.only=T) # load required package
	} # end else statement
} #END packages_used for-loop
if (!all_packages_installed) stop("Need to install required packages")


# "README.md" file documents packages and versions used for future reference
# Function to write lines to the README.md file.  
# Lines have to end with two spaces to cause a carriage return.

write_lines("", paste0("./code/", "README.md"), sep="  \n")

addToREADME <- function(strVec, append=T){
	# Add text with two spaces to get carriage return
  	write_lines(strVec, paste0("./code/", "README.md"), sep="  \n", append=append)
} # end readme function

# addToREADME(paste0("## ", rmarkdown::metadata$title), append=F)

addToREADME(c(date(), ""))
packages_info <- ip[packages_used, c("Package", "Version", "Built")]
addToREADME(c("The packages used in this script are:", "Package, Version, Built"))
apply(packages_info, 1, function(vec) addToREADME(paste(vec, collapse=" ")))
addToREADME("")


# Note define a cluster using makeCluster for FORK parallel backend do not work 

mc.cores <- 20  # Number of cores for parallel processing

# Run scripts to load functions  to simulated gxe data  and fit GS model
source("./code/gxeSim.R")
source("./code/gsModelFit.R") # invoke function to fit GBLUP model
source("./code/varCompModel.R") # to estimate heritability

# scripts to invoke functions to simulate conventional and  GS-based breeding programs
source("./code/convProgram.R") # simulate conventional breeding program
source("./code/broadAdaptation.R")
source("./code/narrowAdaptation.R")

# A function to return  phenotypic value at a p-value of target env.
gvADG <- function(pop, pval){
	return(setPheno(pop, varE=0, p=pval, onlyPheno=T))
} # end gvADG function


# A function to return variance of genetic mean at p-value of target env. 
varG_ADG <- function(pop, pval){
	genoEnv <- setPheno(pop, varE=0, p=pval, onlyPheno=T)
	return(mean(genoEnv^2) - mean(genoEnv)^2)
}


(start.time <- Sys.time())


# Simulate 20-year of burn-in phase for each simulation run or replicate
# long enough to complete a couple breeding cycles so that the parents
# are at least a couple generations removed from the founders

nSimRun <- 50  # Number of simulation repetitions

# executing burnin-phase for a number of simulation replicates (nSimRun)

for(REP in 1:nSimRun){
	cat("\n Simulation run:",REP,"of", nSimRun,"replications", "\n")

  	# Run burn-in phase by advancing years of breeding
  	source("./code/burnin.R") # advance evaluation stages by years
}


# add the simulation parameters to readMe
addToREADME(c(
	paste("The simulation parameters used in this study include:"),
	paste("The number of cores allocated for the simulation job on the server is", mc.cores),
	paste("The number of simulation replicates is", nSimRun),
	paste("The number of breeding cycle in burnin phase is", burninYears),
	paste("The number of breeding cycle in future evaluation phase is", futureYears),
	paste("The number of parental founders for the simulation is", nParents),
 	paste("The number of crosses is", nCrosses),
	paste("The number of progeny per cross is", nProgeny),

  	paste("There are", nQTL, "QTL per chromosome"),
        paste("There are", nSNP, "SNPS chips per chromosome"),
	paste("There are", nQTL+nSNP, "segregating sites per chromosome"),

	paste("Mean oF additive effect is", genMean),
	paste("Variance of additive effect is", varGen),
	paste("Mean of dominance degree is", ddMean),
	paste("Variance of dominance degree  is", ddVar),
	paste("Variance of GxE is", varGxE),"")
 )# end readme parameter list

# Evaluation phase of different breeding scenarios
# mclapply splits these iterations into multiple processes


# This script implemented narrow adaptation breeding program
# with GS to estimate BV in CET based on genomic data to boost
# the selection accuracy at this early stage where phenotype
# information is limited
# GS is used to advance from CET, PYT, AYT, and UYT


# set the seed for reproducible result
set.seed(1234, kind = "L'Ecuyer-CMRG")


#  Scenario 1: conventional breeding program

result.conv <- do.call(rbind,mclapply(1:nSimRun, FUN = baseMod,
	mc.preschedule = TRUE, mc.set.seed = TRUE, 
	mc.cores = mc.cores,mc.silent = FALSE))

print(result.conv)

# The broad  adaptation breeding program
# with GS estimates BV at CET to late stage based on genomic data to boost
# the selection accuracy at this early stage where phenotype
# information is limited
# GS is used to advance from CET, PYT, AYT, and UYT


# Scenario 2: Broad adaptation breeding program

result.broad <- do.call(rbind,mclapply(1:nSimRun, FUN=gsBroad,
	mc.preschedule = TRUE,mc.set.seed = TRUE, 
	mc.cores = mc.cores, mc.silent = FALSE))

print(result.broad)


# Scenario 3: Narrow adaptation - megaEnv1

# This script implemented narrow adaptation breeding program
# with GS to estimate BV in CET based on genomic data to boost
# the selection accuracy at this early stage where phenotype
# information is limited
# GS is used to advance from CET, PYT, AYT, and UYT

result.me1 <- do.call(rbind,mclapply(1:nSimRun,
	FUN = gsNarrow,ME = "ME1",
    	pvalUYT = allLocPvals[c(6,7,8,9)],locUYT= c(6,7,8,9),
	pvalAYT = allLocPvals[c(6,7)], locAYT= c(6,7),
    	pvalPYT = allLocPvals[7], locPYT = 7,
	mc.preschedule = TRUE, mc.set.seed = TRUE, 
	mc.cores = mc.cores, mc.silent = FALSE))

print(result.me1)


# Scenario 3: Narrow adaptation - megaEnv2
# This script implemented narrow adaptation breeding program for ME2
# with GS to estimate BV in CET based on genomic data to boost
# the selection accuracy at this early stage where phenotype
# information is limited
# GS is used to advance individuals from CET, PYT, AYT, and UYT


cat("\n Modeling scenario 3 narrow adaptation - ME2 enabled GS \n")

result.me2 <- do.call(rbind,mclapply(1:nSimRun,
	FUN = gsNarrow,ME = "ME2",
    	pvalUYT = allLocPvals[c(1, 2, 3, 4)], locUYT= c(1,2,3,4),
	pvalAYT = allLocPvals[c(3, 4)],locAYT= c(3,4),		    
	pvalPYT = allLocPvals[3],locPYT = 3,
    	mc.preschedule = TRUE, mc.set.seed = TRUE, 
    	mc.cores = mc.cores,mc.silent = FALSE))


print(result.me2)


######## rowbind, rescale and Saving the simulated data ############

simData <- rbind(result.conv,result.broad,result.me1,result.me2)


# rescale cycle year so that last burnin year is zero

#simData$year <- -(burninYears-1):burninYears # for equal years of burnin and future phase
simData$year <- -(burninYears-1):(burninYears-10) # rescale cycle= 20 burnin + 10 future

# save the simulated data
#write.csv(simData,file=paste0("./data/simData",".csv"),row.names=FALSE)
saveRDS(simData,file=paste0("./data/simData",".rds"))
write.csv(simData,file=paste0("./data/simDataGxE",varGxE,".csv"),row.names=FALSE)


# compute genetic parameters from selection cycles of the simulated data
simSumm <- simData %>%
	dplyr::select(-c(nParPYT,nParAYT,nParUYT)) %>%
  	dplyr::filter(year >= 0) %>%
  	group_by(scenario,year) %>%
  	dplyr::summarise(
    	genMean = mean(meanGV),
    	seGenMean = sd(meanGV)/sqrt(length(simRun)),
    	genVar = mean(varGV),
    	seGenVar = sd(varGV)/sqrt(length(simRun)),
	selAccur= mean(accPYT),
    	seAccur = sd(accPYT)/sqrt(length(simRun)),.groups = 'drop')


simSumm$scenario <- factor(simSumm$scenario,
	levels = c("Conv","Broad adaptation","Narrow adaptation"),
        labels = c("Conventional","Broad adaptation","Narrow adaptation"))


write.csv(simSumm,file=paste0("./data/simSumm",varGxE,".csv"),row.names=FALSE)


# fit linear regression of genetic gain on the breeding cycle of evaluation phase
# to calculate rate of genetic gain per year

fitConv <- lm(genMean ~ year, data = simSumm[simSumm$scenario=='Conventional' & simSumm$year > 0,])
cat("Rate of genetic gain per year for conventional breeding program is",
	round(coef(fitConv)[2],2),"\n")


fitBroad <- lm(genMean ~ year, data = simSumm[simSumm$scenario=='Broad adaptation' & simSumm$year > 0,])
cat("Rate of genetic gain per year for broad adaptation breeding program is",
	round(coef(fitBroad)[2],2),"\n")


fitNarrow <- lm(genMean ~ year, data = simSumm[simSumm$scenario=='Narrow adaptation' & simSumm$year > 0,])
cat("Rate of genetic gain per year for narrow adaptation breeding program is",
    	round(coef(fitNarrow)[2],2),"\n")

(end.time <- Sys.time())
(end.time - start.time)

# clean the folder by deleting all the burnin files (*.rda)
unlink("./data/bu*.rda")


# plot distribution of the number of individuals that constitute the parental candidate

parCount <- simData %>%
        dplyr::filter(year > 0) %>%
        dplyr::select(simRun,year,scenario,nParPYT,nParAYT,nParUYT) %>%
	gather(key = stage,value = nPar,-c(simRun,year,scenario))

parCount$scenario <- factor(parCount$scenario,
	levels = c("Conv","Broad adaptation","Narrow adaptation"),
	labels=c("Conventional","Broad adaptation","Narrow adaptation"))


parCount$stage <- factor(parCount$stage,
	levels = c("nParPYT","nParAYT","nParUYT"),
        labels=c("PYT","AYT","UYT"))

distParCand  <- ggplot(data = parCount, aes(x=scenario,y=nPar, fill=stage)) +
	geom_boxplot( stat = "boxplot",  outlier.colour = "red", outlier.shape=16,
        outlier.size = 0.5, na.rm=TRUE) +
  	labs(x= "Breeding strategy", y= "Parental lines") + theme_bw() +
  	theme(axis.title = element_text(colour="black", size=12),
        	plot.title = element_text(hjust = 0.5,lineheight=.5,colour="black", size=12),
         	axis.text = element_text(face="bold", size=7))

# print(distParCand)
# save the plot to a file
ggsave("./output/distParentalCandidate.jpeg",height=4.5, width=6.5, units="in", dpi=300)


# A scatter plot of rate of average genetic value and cycle-year with regression equation line

rateGenGain <- ggplot(data=simSumm[simSumm$year>0,], aes(x = year, y = genMean,color=scenario)) +
	geom_point() + labs(x="Cycles of selection",y="Genetic mean")+
  	stat_smooth(aes(fill = scenario, color = scenario), method = "lm", se=FALSE,
      	formula = y ~ poly(x, 1, raw = TRUE)) +
  	stat_regline_equation(
    	aes(label =  paste(after_stat(eq.label), after_stat(adj.rr.label), sep = "~~~~")),
    	formula = y ~ poly(x, 1, raw = TRUE)) + theme_bw()

# save the plot to a file
ggsave("./output/rateGeneticGain.jpeg",height=4.2, width=6.5, units="in", dpi=300)



# plot trend of genetic mean

genMeanPlot <- ggplot(data = simSumm,aes(x = year,y = genMean, color = scenario))+
	geom_point()+
  	geom_line(linewidth = 1)+
  	geom_errorbar(aes(ymin = genMean - seGenMean, ymax = genMean + seGenMean), width = 0.2,
	position = position_dodge(0.05))+
  	guides(scale = "none")+
  	theme_bw()+ theme(legend.position = c(0.03,0.98),
        legend.justification = c("left","top"))+
  	scale_x_continuous("Cycles of selection",limits = c(0,futureYears))+
  	scale_y_continuous("Genetic mean",limits = c(0,NA))

# print(genMeanPlot)

# save the plot to a file
ggsave("./output/geneticMean.jpeg",height = 4.5, width = 6.5, units = "in", dpi = 300)



# plot trend of genetic variance
genVarPlot <- ggplot(data = simSumm,aes(x = year,y = genVar,color = scenario))+ geom_point()+
	geom_line(linewidth = 1)+
  	geom_errorbar(aes(ymin = genVar-seGenVar, ymax=genVar+seGenVar), width = 0.2,
                position = position_dodge(0.05))+
  	guides(scale = "none")+
  	theme_bw()+ theme(legend.position = c(0.03,0.05),
        legend.justification = c("left","bottom"))+
  	scale_x_continuous("Cycles of selection",limits=c(0,futureYears))+
  	scale_y_continuous("Genetic variance",limits=c(0,NA))

# print(genVarPlot)

# save the plot to a file
ggsave("./output/geneticVariance.jpeg",height = 4.5, width = 6.5, units = "in", dpi = 300)



# boxplot of comparing the distribution of  accuracy of selection criterion acros the stages

selAccur <- simData %>%
	dplyr::filter(year>0) %>%
	dplyr::select(simRun, year,scenario,accSDN,accCET,accPYT,accAYT,accUYT) %>%
  	gather(key = stage,value = accur,-c(simRun,year,scenario))

selAccur$scenario <- factor(selAccur$scenario,
	levels = c("Conv","Broad adaptation","Narrow adaptation"),
	labels=c("Conventional","Broad adaptation", "Narrow adaptation"))

selAccur$stage <- factor(selAccur$stage,
        levels = c("accSDN","accCET","accPYT","accAYT","accUYT"),
        labels = c("SDN","CET","PYT","AYT","UYT"))

distSelAccur <- ggplot(data = selAccur, aes(x = scenario,y = accur, fill = stage)) +
	geom_boxplot( stat = "boxplot",  outlier.colour = "red", outlier.shape = 16,
        outlier.size = 0.5, na.rm=TRUE) +
  	labs(x= "Breeding strategy", y= "Selection accuracy") + theme_bw() +
  	theme(axis.title = element_text(colour="black", size=12),
   	plot.title = element_text(hjust = 0.5,lineheight = .5,colour="black", size = 12),
  	axis.text = element_text(face = "bold", size = 7))

# print(distSelAccur)
# save the plot to a file
ggsave("./output/distSelectionAccuracy.jpeg",height = 4.5, width = 6.5, units = "in", dpi = 300)


# plot Mean of selection accuracy
selAccur <- simData %>%
	dplyr::filter(year>0) %>%
  	group_by(scenario) %>%
  	dplyr::summarise(
		accSDN=mean(accSDN),
    		accCET=mean(accCET),
    		accPYT=mean(accPYT),
    		accAYT=mean(accAYT),
    		accUYT=mean(accUYT)) %>%
		gather(key = stage,value = accur,-c(scenario))

	selAccur$scenario <- factor(selAccur$scenario,
          	levels = c("Conv","Broad adaptation","Narrow adaptation"),
                labels = c("Conventional","Broad adaptation","Narrow adaptation"))

	selAccur$stage <- factor(selAccur$stage,
                levels = c("accSDN","accCET","accPYT","accAYT","accUYT"),
                labels=c("SDN","CET","PYT","AYT","UYT"))

# print(selAccur)

selAccPlot <- ggplot(data = selAccur,aes(x = scenario,y = accur,fill = stage)) +
geom_bar(stat = "identity", position = "dodge")+
  geom_text(aes(x = scenario,y = accur,label = round(accur,2)),
            position = position_dodge(width = 0.9),vjust=-0.25)+
  labs(x = "Breeding programme", y = "Selection accuracy")

#print(selAccPlot)

 # save the plot to a file
ggsave("./output/selAccuracyMean.jpeg",height=4.2, width=6.5, units="in", dpi=300)


# plot line trends of the accuracy at PYT stage

selAccurPlot <- ggplot(data=simSumm,aes(x=year,y=selAccur,color=scenario))+
  geom_point()+
  geom_line(linewidth = 1)+
  geom_errorbar(aes(ymin=selAccur-seAccur, ymax=selAccur+seAccur), width=0.2,
                position=position_dodge(0.05))+
  guides(scale="none")+
  theme_bw()+
  theme(legend.position = c(0.95,0.10),
        legend.justification = c("right","bottom"))+
  scale_x_continuous("Cycles of selection",limits=c(0,futureYears))+
  scale_y_continuous("Selection accuracy",limits=c(0,NA))

# print(selAccurPlot)

# save the plot to a file
ggsave("./output/selectionAccuracyTrend.jpeg",height=4.2, width=6.5, units="in", dpi=300)



# plot both  narrow and broad-sense  heritability
heritab <- simData %>%
  dplyr::filter(year>0) %>%
  group_by(scenario) %>%
  dplyr::summarise(h2=mean(h2, na.rm=TRUE),H2 =mean(H2),.groups = 'drop') %>%
  gather(key = heritability,value = estimate,-c(scenario))

heritab$scenario <- factor(heritab$scenario,
                          levels = c("Conv","Broad adaptation","Narrow adaptation"),
                          labels = c("Conventional","Broad adaptation","Narrow adaptation"))


heritabPlot <- ggplot(data = heritab,aes(x = scenario,y = estimate,fill = heritability)) +
  geom_bar(stat = "identity", position = "dodge")+
  geom_text(aes(x = scenario,y = estimate,label = round(estimate,2)),
            position = position_dodge(width = 0.9),vjust=-0.25)+
  labs(x = "Breeding strategy", y = "Heritability estimate") + theme_bw() +
  scale_fill_discrete(name = "Heritability", 
                      labels = c("Narrow-sense", "Broad-sense"))


#print(heritabPlot)

 # save the plot to a file
ggsave("./output/heritabEstimate.jpeg",height=4.5, width=6.5, units="in", dpi=300)

broadH2 <- heritab %>%
  dplyr::filter(heritability=="H2")

H2Plot <- ggplot(data = broadH2,aes(x = scenario,y = estimate,fill = scenario)) +
        geom_bar(stat = "identity", position = "dodge")+
        geom_text(aes(x = scenario,y = estimate,label = round(estimate,2)),
        position = position_dodge(width = 0.9),vjust=-0.25)+
        labs(x = "Breeding strategy", y = "Heritability estimate")


#print(H2Plot)
# save the plot to a file
ggsave("./output/broad_SenseH2Summary.jpeg",height = 4.5, width = 6.5, units = "in", dpi = 300)

addToREADME(c(
        paste("The total time for running the simulation is ",(end.time - start.time))
))
