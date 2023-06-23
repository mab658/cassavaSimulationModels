# Plot the result
rm(list=ls())
library(ggplot2)
library(tidyr)
library(dplyr)
library(ggpubr)

# working folder
setwd("/Users/mab658/Documents/cassavaSimulationModels/simulationPlot")

simDataGxE0 <- read.csv(file="simDataGxE0.csv",header = T)
simDataGxE2 <- read.csv(file="simDataGxE4.csv",header = T)

simData <- rbind(simDataGxE0,simDataGxE4)

# plot distribution of the number of individuals that constitute the parental candidate

simData$varGxE <- factor(simData$varGxE, 
                          labels = c("varGxE = 0.0", "varGxE = 2.0"))

parCount <- simData %>%
  dplyr::filter(year > 0 ) %>%
  dplyr::select(simRun,year,scenario, varGxE,nParPYT,nParAYT,nParUYT) %>%
  gather(key = stage,value = nPar,-c(simRun,year,scenario, varGxE))

parCount$scenario <- factor(parCount$scenario,
                            levels = c("Conv","Broad adaptation",
                                       "Narrow adaptation"),
                            labels=c("Conventional","Broad adaptation",
                                     "Narrow adaptation"))

parCount$stage <- factor(parCount$stage,
                         levels = c("nParPYT","nParAYT","nParUYT"),
                         labels=c("PYT","AYT","UYT"))


distParCand  <- ggplot(data = parCount, aes(x=scenario,y=nPar, fill=stage)) +
  geom_boxplot( stat = "boxplot",  outlier.colour = "red", outlier.shape=16,
                outlier.size = 0.5, na.rm=TRUE) +
  labs(x= "Breeding strategy", y= "# of parental candidate") + theme_bw() +
  theme(axis.title = element_text(colour="black", size=12),
        plot.title = element_text(hjust = 0.5,lineheight=.5,colour="black", size=12),
        axis.text = element_text(face="bold", size=7)) +
  theme(legend.position = "bottom",
        legend.justification = "center")+
  facet_grid(.~varGxE,scales = "free")

 print(distParCand)
# save the plot to a file
ggsave("./output/distParentalCandidate.jpeg",height=4.2, width=6.5, units="in", dpi=300)


# compute genetic parameters from selection cycles of the simulated data

simSumm <- simData %>%
  dplyr::select(-c(nParPYT,nParAYT,nParUYT)) %>%
  dplyr::filter(year >= 0) %>%
  group_by(scenario, varGxE,year) %>%
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


# A scatter plot of average genetic value and cycle-year with regression equation line

rateGenGain <- ggplot(data=simSumm[simSumm$year>0,], aes(x = year, y = genMean,color=scenario)) +
  geom_point() +labs(y="Genetic mean")+
  stat_smooth(aes(fill = scenario, color = scenario), method = "lm", se=FALSE,
              formula = y ~ poly(x, 1, raw = TRUE)) +
  stat_regline_equation(
    aes(label =  paste(after_stat(eq.label), after_stat(adj.rr.label), sep = "~~~~")),
    formula = y ~ poly(x, 1, raw = TRUE)
  ) + labs(x="Cycles of selection") + theme_bw() +
  theme(legend.position = "bottom",
    legend.justification = "center")+
  scale_x_continuous("Cycles of selection",limits = c(0,10))+ #10=futuYears
  scale_y_continuous("Genetic mean",limits = c(0,NA))+
  
  facet_grid(.~ varGxE,scales = "free")
  
 
print(rateGenGain) 


# save the plot to a file
ggsave("./output/rateGeneticGain.jpeg",height=4.2, width=6.5, units="in", dpi=300)


# plot trend of genetic mean

genMeanPlot <- ggplot(data = simSumm,aes(x = year,y = genMean, color = scenario))+
  geom_point()+
  geom_line(linewidth = 1)+
  geom_errorbar(aes(ymin = genMean - seGenMean, ymax = genMean + seGenMean), width = 0.2,
                position = position_dodge(0.05))+
  guides(scale = "none")+
  theme_bw()+
  theme(legend.position = c(0.03,0.98),
        legend.justification = c("left","top"))+
  scale_x_continuous("Cycles of selection",limits = c(0,10))+ #10=futuYears
  scale_y_continuous("Genetic mean",limits = c(7,NA))+
  theme(legend.position = "bottom",
        legend.justification = "center")+
  facet_wrap(.~varGxE,scales = "free")

print(genMeanPlot)

# save the plot to a file
ggsave("./output/geneticMean.jpeg",height = 4.2, width = 6.5, units = "in", dpi = 300)




# plot trend of genetic variance
genVarPlot <- ggplot(data = simSumm,aes(x = year,y = genVar,color = scenario))+ geom_point()+
  geom_line(linewidth = 1)+
  geom_errorbar(aes(ymin = genVar-seGenVar, ymax=genVar+seGenVar), width = 0.2,
                position = position_dodge(0.05))+
  guides(scale = "none")+
  theme_bw()+
  theme(legend.position = "bottom",
        legend.justification = "center")+
  scale_x_continuous("Cycles of selection",limits=c(0,10))+ #NA=futureYear
  scale_y_continuous("Genetic variance",limits=c(0,NA)) +
  theme(legend.position = "bottom",
        legend.justification = "center")+
  facet_wrap(.~varGxE,scales = "free")

print(genVarPlot)

# save the plot to a file
ggsave("./output/geneticVariance.jpeg",height = 4.2, width = 6.5, units = "in", dpi = 300)


# boxplot of comparing the distribution of  accuracy of selection criterion acros the stages

selAccur <- simData %>%
  dplyr::filter(year>0) %>%
  dplyr::select(simRun, year,scenario,varGxE,accSDN,accCET,accPYT,accAYT,accUYT) %>%
  gather(key = stage,value = accur,-c(simRun,year,scenario,varGxE))

selAccur$scenario <- factor(selAccur$scenario,
                            levels = c("Conv","Broad adaptation",
                                       "Narrow adaptation"),
                            labels=c("Conventional","Broad adaptation",
                                     "Narrow adaptation"))

selAccur$stage <- factor(selAccur$stage,
                         levels = c("accSDN","accCET","accPYT","accAYT","accUYT"),
                         labels = c("SDN","CET","PYT","AYT","UYT"))

distSelAccur <- ggplot(data = selAccur, aes(x = scenario,y = accur, fill = stage)) +
  geom_boxplot( stat = "boxplot",  outlier.colour = "red", outlier.shape = 16,
                outlier.size = 0.5, na.rm=TRUE) +
  labs(x= "Breeding strategy", y= "Selection accuracy") + theme_bw() +
  theme(axis.title = element_text(colour="black", size=12),
        plot.title = element_text(hjust = 0.5,lineheight = .5,colour="black", size = 12),
        axis.text = element_text(face = "bold", size = 7))+
  theme(legend.position = "bottom",
        legend.justification = "center")+
  facet_grid(.~varGxE,scales = "free")

print(distSelAccur)
# save the plot to a file
ggsave("./output/distSelectionAccuracy.jpeg",height=4.2, width=6.5, units="in", dpi=300)



# plot Mean of selection accuracy
selAccur <- simData %>%
  dplyr::filter(year>0) %>%
  group_by(scenario, varGxE) %>%
  dplyr::summarise(
    accSDN=mean(accSDN),
    accCET=mean(accCET),
    accPYT=mean(accPYT),
    accAYT=mean(accAYT),
    accUYT=mean(accUYT),.groups = 'drop') %>%
  gather(key = stage,value = accur,-c(scenario, varGxE))

selAccur$scenario <- factor(selAccur$scenario,
                            levels = c("Conv","Broad adaptation","Narrow adaptation"),
                            labels = c("Conventional","Broad adaptation","Narrow adaptation"))

selAccur$stage <- factor(selAccur$stage,
                         levels = c("accSDN","accCET","accPYT","accAYT","accUYT"),
                         labels=c("SDN","CET","PYT","AYT","UYT"))

selAccPlot <- ggplot(data = selAccur,aes(x = scenario,y = accur,fill = stage)) +
  geom_bar(stat = "identity", position = "dodge")+
  geom_text(aes(x = scenario,y = accur,label = round(accur,2)),
            position = position_dodge(width = 0.9),vjust=-0.25,size=2.5)+
  labs(x = "Breeding strategy", y = "Selection accuracy") + theme_bw() +
  theme(axis.title = element_text(colour="black", size=12),
        plot.title = element_text(hjust = 0.5,lineheight=.5,colour="black", size=12),
        axis.text = element_text(face="bold", size=7)) +
  theme(legend.position = "bottom",
        legend.justification = "center")+
  facet_grid(.~varGxE,scales = "free")

#print(selAccPlot)

# save the plot to a file
ggsave("./output/meanSelectionAccuracy.jpeg",height=4.2, width=6.5, units="in", dpi=300)



# plot both narrow and broad-sense heritability

heritab <- simData %>%
  dplyr::filter(year>0) %>%
  group_by(scenario,varGxE) %>%
  dplyr::summarise(h2=mean(h2, na.rm=TRUE),H2= mean(H2),.groups = 'drop') %>%
  gather(key = heritability,value = estimate,-c(scenario,varGxE))

heritab$scenario <- factor(heritab$scenario,
                          levels = c("Conv","Broad adaptation","Narrow adaptation"),
                          labels = c("Conventional","Broad adaptation","Narrow adaptation"))


heritPlot <- ggplot(data = heritab,aes(x = scenario,y = estimate,fill = heritability)) +
  geom_bar(stat = "identity", position = "dodge")+
  geom_text(aes(x = scenario,y = estimate,label = round(estimate,2)),
            position = position_dodge(width = 0.9),vjust=-0.25)+
  labs(x = "Breeding strategy", y = "Heritability estimate") + theme_bw() +
  scale_fill_discrete(name = "Heritability", 
                      labels = c("Narrow-sense", "Broad-sense"))+
  theme(axis.title = element_text(colour="black", size=12),
        plot.title = element_text(hjust = 0.5,lineheight = .5,colour="black", size = 12),
        axis.text = element_text(face = "bold", size = 7))+
  theme(legend.position = "bottom",
        legend.justification = "center")+
  facet_grid(.~varGxE,scales = "free")


plot(heritPlot)
ggsave("./output/heritabEstimate.jpeg",height=4.5, width=6.5, units="in", dpi=300)



# Plot of only broad-sense heritability

broadH2 <- heritab %>%
  dplyr::filter(heritability=="H2")

H2Plot <- ggplot(data = broadH2,aes(x = scenario,y = estimate,fill = scenario)) +
  geom_bar(stat = "identity", position = "dodge")+
  geom_text(aes(x = scenario,y = estimate,label = round(estimate,2)),
            position = position_dodge(width = 0.9),vjust=-0.25)+
  labs(x = "Breeding strategy", y = "Broad-sense heritability") + theme_bw()+
  theme(axis.title = element_text(colour="black", size=12),
        plot.title = element_text(hjust = 0.5,lineheight = .5,colour="black", size = 12),
        axis.text = element_text(face = "bold", size = 7),
        legend.position = "bottom",
        legend.justification = "center")+
  facet_grid(.~varGxE,scales = "free")

print(H2Plot)
# save the plot to a file
ggsave("./output/broad_Sense_heritab.jpeg",height = 4.5, width = 6.5, units = "in", dpi = 300)


# plot line trends of the accuracy at PYT stage

selAccurPlot <- ggplot(data=simSumm,aes(x=year,y=selAccur,color=scenario))+
  geom_point()+
  geom_line(size=1)+
  geom_errorbar(aes(ymin=selAccur-seAccur, ymax=selAccur+seAccur), width=0.2,
                position=position_dodge(0.05))+
  guides(scale="none")+
  theme_bw()+
  theme(legend.position = "bottom",
        legend.justification = "center")+
  scale_x_continuous("Year",limits=c(0,10))+
  scale_y_continuous("Selection accuracy",limits=c(0,NA))+
  facet_grid(.~varGxE,scales = "free")

# print(selAccurPlot)

# save the plot to a file
ggsave("./output/PYTselAccuracyTrend.jpeg",height=4.2, width=6.5, units="in", dpi=300)


# Comparing distribution of additive and non-additive effect
genEff <- simData %>%
  dplyr::filter(year > 0 & scenario != "Conv") %>%
  dplyr::select(scenario,varGxE,ebv,egv) %>%
  gather(key = effect,value = estimate,-c(scenario, varGxE))

distrEff <- ggplot(genEff, aes(x=estimate, colour=effect, fill=effect)) + geom_density()+
  labs(x="Genetic effect")+ theme_bw()+
  scale_fill_discrete(name = "Effect",
                      labels = c("Additive effect", "Non-additive effect")) +
  facet_grid(.~varGxE,scales = "free")
  
ggsave("./output/distrGenEffect.jpeg",height = 4.5, width = 6.5, units = "in", dpi = 300)






