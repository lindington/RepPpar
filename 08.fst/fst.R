###Populationwise Fst

library(ade4)
library(ggplot2)
library(extrafont)
library(RColorBrewer)
library(reshape2)
library(tidyverse)
library(RColorBrewer)

#Reading in files and creating the matrices
getwd()


por_PopwiseFst <- read.csv(file = "por_PopulationwiseFst.csv", header = T, row.names = 1, sep = ",", dec = ".")
por_PopwiseDist <- read.csv(file = "porPopulationwiseDistance.csv", header = T, row.names = 1, sep = ",", dec = ".")

por_PopFstMatrix <- dist(por_PopwiseFst)
por_PopDistMatrix <- dist(por_PopwiseDist)

as.matrix(por_PopFstMatrix)
as.matrix(por_PopDistMatrix)

#Performing Mantel's test

mantel.rtest(por_PopFstMatrix, por_PopDistMatrix, nrepet = 1000)

#Scatterplot

por_PopFstXDistDF <- read.csv(file = "por_fstxdist.csv", header = F, row.names = 1, dec = ".", col.names = c("Pair","Fst","Distance.km"))


#Gray area is 0.95 confidence interval

ggplot(data = por_PopFstXDistDF, aes(x = Distance.km, y = Fst)) +
  geom_point() +
  geom_smooth(method='lm') +
  expand_limits(y=0) +
  scale_y_continuous(expand = c(0.01, 0.01), limits = c(0, 0.54)) +
  scale_x_continuous(expand = c(0.01, 0.01), limits = c(0.0, 21.1)) +
  ggtitle("Portalet Transect: Isolation by distance") +
  theme(plot.title = element_text(size = 24)) +
  ylab("Fst") +
  xlab("Distance") +
  theme(text=element_text(size=16))+
  theme(axis.text.x = element_text(size = 8, vjust = 0.5, hjust=0.45))+
  theme(axis.text.y = element_text(size = 8, vjust = 0.5, hjust=0.45))


###Fst France

FranceFst <- PopwiseFst[1:6,1:6]
FranceDist <- PopwiseDist[1:6,1:6]

FranceFstMatrix <- dist(FranceFst)
FranceDistMatrix <- dist(FranceDist)

as.matrix(FranceFstMatrix)
as.matrix(FranceDistMatrix)

#Performing Mantel's test

mantel.rtest(FranceFstMatrix, FranceDistMatrix, nrepet = 1000)

#Scatterplot

FranceFstXDistDF <- read.csv(file = "fstxdist_gab.txt", header = F, row.names = 1, sep = ",", dec = ".", col.names = c("Pair","Fst","Distance.km"))

#Gray area is 0.95 confidence interval

ggplot(data = FranceFstXDistDF, aes(x = Distance.km, y = Fst)) +
  geom_point() +
  geom_smooth(method='lm') +
  expand_limits(y=0) +
  scale_y_continuous(expand = c(0.01, 0.01), limits = c(-0.006, 0.325)) +
  scale_x_continuous(expand = c(0.01, 0.01), limits = c(0.0, 10.6)) +
  ggtitle("Isolation by distance north of hybridzone") +
  theme(plot.title = element_text(size = 24)) +
  ylab("Fst") +
  xlab("Distance") +
  theme(text=element_text(size=16))+
  theme(axis.text.x = element_text(size = 8, vjust = 0.5, hjust=0.45))+
  theme(axis.text.y = element_text(size = 8, vjust = 0.5, hjust=0.45))


###Fst Spain

SpainFst <- PopwiseFst[8:11,8:11]
SpainDist <- PopwiseDist[8:11,8:11]

SpainFstMatrix <- dist(SpainFst)
SpainDistMatrix <- dist(SpainDist)

as.matrix(SpainFstMatrix)
as.matrix(SpainDistMatrix)

#Performing Mantel's test

mantel.rtest(SpainFstMatrix, SpainDistMatrix, nrepet = 1000)

#Scatterplot

SpainFstXDistDF <- read.csv(file = "fstxdist_esc.txt", header = F, row.names = 1, sep = ",", dec = ".",col.names = c("Pair","Fst","Distance.km"))

#Gray area is 0.95 confidence interval

ggplot(data = SpainFstXDistDF, aes(x = Distance.km, y = Fst)) +
  geom_point() +
  geom_smooth(method='lm') +
  expand_limits(y=0) +
  scale_y_continuous(expand = c(0.01, 0.01), limits = c(-0.008, 0.125)) +
  scale_x_continuous(expand = c(0.01, 0.01), limits = c(0.0, 3.4)) +
  ggtitle("Isolation by distance south of hybridzone") +
  theme(plot.title = element_text(size = 24)) +
  ylab("Fst") +
  xlab("Distance") +
  theme(text=element_text(size=16))+
  theme(axis.text.x = element_text(size = 8, vjust = 0.5, hjust=0.45))+
  theme(axis.text.y = element_text(size = 8, vjust = 0.5, hjust=0.45))
