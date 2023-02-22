#############################
# Individual Heterozygosity #
#############################

#Importing packages and Avenir fonts
library(methods)
library(ggplot2)
library(scales)
library(reshape2)
library(ggpubr)
library(dplyr)
library(viridis)

setwd("C:/Users/Jag/Documents/LMU/work/01.capture/05.hz")
#Reading in the data
IndHZData <- read.csv(file = "219_sum_sfs_ind.ml", header = F, sep = "", dec = ".")
IndHZData$pop <- c(rep("OUT",10),rep("GAB",5),rep("HER",10),rep("SOQ",10),rep("TOU",10),rep("ARA",10),rep("POR",13),rep("MUL",7),rep("SOQ",1),rep("MUL",6),rep("TRO",13),rep("FOR",10),rep("PAZ",10),rep("LAN",10),rep("ESC",5),rep("STP",10),rep("AIN",10),rep("URD",10),rep("OTS",10),rep("ARI",10),rep("LEK",9),rep("MAR",10),rep("ALM",10),rep("VEN",10)) 
IndHZData$colrs <- c(rep("black",10),rep("#A6CEED",45),rep("#FAD9BE",19),rep("#A6CEED",1),rep("#FAD9BE",6),rep("#ED797C",49),rep("#ED797C",20),rep("#FAD9BE",20),rep("#ED797C",29),rep("#FAD9BE",20)) 
#Formatting data
IndHZs <- IndHZData[1]
IndHZs[2] <- IndHZData$pop
IndHZs[3] <- IndHZData[3]/sum(IndHZData[2:4])

#Reading in hexadecimal colour vectors
ViridisColours <- c("#18D6CBFF", "#455BCDFF", "#A2FC3CFF", "#F05B12FF", "#30123BFF", "#C42503FF", "#FEA632FF", "#7A0403FF", "#E1DD37FF", "#46F884FF", "#3E9BFEFF")

#Reading in order data
LocalityOrder <- factor(IndHZs$V2 , levels = c("OUT","GAB", "HER", "SOQ", "TOU", "ARA", "POR", "MUL", "TRO","FOR", "PAZ", "LAN", "ESC", "STP", "AIN", "URD","OTS","ARI","LEK","MAR","ALM","VEN"))
LocalityOrder
LocalityOrder2 <- c("OUT","GAB", "HER", "SOQ", "TOU", "ARA", "POR", "MUL", "TRO","FOR", "PAZ", "LAN", "ESC", "STP", "AIN", "URD","OTS","ARI","LEK","MAR","ALM","VEN")
colourorder<-c(rep("black",1),rep("#a6ceed",3),rep("#FAD9BE",5),rep("#ED797C",4),rep("#a6ceed",2),rep("#FAD9BE",2),rep("#a6ceed",3),rep("#FAD9BE",2))

colperpop<-c("OUT"="black","GAB"="#A6CEED","HER"="#A6CEED","SOQ"="#A6CEED","TOU"="#A6CEED","ARA"="#A6CEED","POR" = "#FAD9BE", "MUL"="#FAD9BE", "TRO"="#FAD9BE","FOR"="#ED797C","PAZ"="#ED797C", "LAN"="#ED797C", "ESC"="#ED797C", "STP"="#A6CEED", "AIN"="#A6CEED", "URD"= "#FAD9BE","OTS"= "#FAD9BE","ARI"="#A6CEED","LEK"="#A6CEED","MAR"="#A6CEED","ALM"= "#FAD9BE","VEN"="#FAD9BE")
###Per Individual Heterozygosity boxplot
IndHZs$colrs<-as.character(IndHZData$colrs)

ggplot(IndHZs, aes(x = factor(V2, levels=LocalityOrder2),V3)) +
  geom_boxplot(aes(fill = LocalityOrder), outlier.shape = NA) +
  geom_jitter(width = 0.075, height = 0, size = 1.96) +
  ggtitle("Individual Heterozygosity Per Population") +
  theme(plot.title = element_text(size = 18)) +
  ylab("Heterozygosity") +
  xlab("Locality") +
  theme(text=element_text(size=16)) +
  theme(axis.text.x = element_text(size = 10)) + 
  scale_fill_manual(values=colourorder, name="Locality",breaks=LocalityOrder2,labels=c("Italy(out)","Gabas","Hermine","Soques","Tourmont","Araille","Portalet","Mulas","La Troya","Formigal","Pazino","Lanuza","Escarilla","St Pee","Ainoha","Urdax","Otsondo","Arizcun","Lekaroz","Mardea","Almandoz","Ventas de Arraiz"))

###Min, max and avg heterozygosity
which.min(IndHZs[,3])
which.max(IndHZs[,3])
mean(IndHZs[,3])
