############################################################################################
###################### PLOT THE CLINES FROM THE PREPROCESSED DATASETS ######################
############################################################################################
### load packages

install.packages("scico")
library(tidyverse)
library(ggplot2)
library(data.table)
library(hzar)
library(dplyr)
library(scico)


### set to working directory
getwd()
setwd("C:/Users/Jag/Documents/LMU/work/01.capture/07.hzar")

## read in summary file 
all_frqs.summary <- readRDS("all_frqs_summary")


## get snp names

all_frqs <- read.csv("../df_frqs.csv")

snplist<- seq(3,ncol(all_frqs),2)
loclist<-c(colnames(all_frqs)[snplist])

## get delta freq change

# polarise: ascending

for (i in snplist){
  if (all_frqs[1,i] > all_frqs[11,i]){
    all_frqs[[i]]<- sapply(all_frqs[[i]], function(x) 1-x)
  }
}

deltaf<-list()
for (i in snplist){
  k<-((i-1)/2)
  deltaf[[k]]<-all_frqs[11,i]-all_frqs[1,i]
}

# extract centre, width, credibility interval values
all_frqs.centres<-list()
for (i in 1:length(all_frqs.summary)){
  all_frqs.centres[[i]]<-all_frqs.summary[[i]]$statistics[[1]]
}
all_frqs.width<-list()
for (i in 1:length(all_frqs.summary)){
  all_frqs.width[[i]]<-all_frqs.summary[[i]]$statistics[[2]]
}

all_frqs.centresCIlo<-list()
for (i in 1:length(all_frqs.summary)){
  all_frqs.centresCIlo[[i]]<-all_frqs.summary[[i]]$quantiles[1,1]
}
all_frqs.centresCIhi<-list()
for (i in 1:length(all_frqs.summary)){
  all_frqs.centresCIhi[[i]]<-all_frqs.summary[[i]]$quantiles[1,5]
}
all_frqs.widthCIlo<-list()
for (i in 1:length(all_frqs.summary)){
  all_frqs.widthCIlo[[i]]<-all_frqs.summary[[i]]$quantiles[2,1]
}
all_frqs.widthCIhi<-list()
for (i in 1:length(all_frqs.summary)){
  all_frqs.widthCIhi[[i]]<-all_frqs.summary[[i]]$quantiles[2,5]
}

allsum<-as.data.frame(loclist)
allsum$centre<-unlist(all_frqs.centres)
allsum$centreCIlo<-unlist(all_frqs.centresCIlo)
allsum$centreCIhi<-unlist(all_frqs.centresCIhi)
allsum$width<-unlist(all_frqs.width)
allsum$widthCIlo<-unlist(all_frqs.widthCIlo)
allsum$widthCIhi<-unlist(all_frqs.widthCIhi)
allsum$deltaf<-unlist(deltaf)
allsum$slope<-allsum$deltaf/(allsum$width/3000)

## plot 

ggplot(allsum, aes(centre))+ 
  geom_histogram(binwidth=500, fill="#69b3a2") +
  theme_bw()

ggplot(allsum, aes(width))+ 
  geom_histogram(binwidth=500, fill="#69b3a2") +
  theme_bw()

ggplot(allsum, aes(centre,width))+ 
  geom_point() +
  theme_bw()

ggplot(allsum, aes(centre,centreCI))+ 
  geom_point() +
  theme_bw()

ggplot(allsum, aes(centre,deltaf))+
  geom_point()+
  theme_bw()

### restrict deltaf, finalise table

allsum<-allsum %>%
  filter(deltaf > 0.5)

ggplot(allsum,aes(centre))+
  geom_histogram(binwidth=1000, fill="#69b3a2")+
  theme_bw()

ggplot(allsum,aes(width))+
  geom_histogram(binwidth=1000, fill="#69b3a2")+
  theme_bw()


ggplot(allsum, aes(centre,slope))+ 
  #stat_density2d(aes(fill=..level..), geom="polygon") +
  stat_density2d(aes(fill = ..density..), geom = "raster", contour = FALSE) +
  #scale_fill_gradient(low="#ffc9c9", high="#b00000") +
  scale_fill_viridis_c(option="magma",direction=-1) +
  geom_point(fill="seagreen",alpha=0.8,size=0.3) +
  theme_bw()

ggplot(allsum, aes(centre,width))+ 
  #stat_density2d(aes(fill=..level..), geom="polygon") +
  stat_density2d(aes(fill = ..density..), geom = "raster", contour = FALSE) +
  #scale_fill_gradient(low="#ffc9c9", high="#b00000") +
  scale_fill_viridis_c(option="magma",direction=-1) +
  geom_point(fill="seagreen",alpha=0.8,size=0.3) +
  theme_bw()

allsum <- allsum %>%
  unite('width_CI',widthCIlo:widthCIhi,sep=" - ", remove = TRUE) %>%
  unite('centre_CI',centreCIlo:centreCIhi,sep=" - ", remove = TRUE) %>%
  separate(loclist,into = c("X","SNP"),sep="X", remove = TRUE) 

allsum <- allsum[2:length(allsum)]

allsum$width_CI <- sub("^","[",allsum$width_CI)
allsum$width_CI <- sub("$","]",allsum$width_CI)
allsum$centre_CI <- sub("^","[",allsum$centre_CI)
allsum$centre_CI <- sub("$","]",allsum$centre_CI)

write.csv(allsum,"allsum.csv", row.names = FALSE)

## find number of snps outside width (based on HI clines):

allsum$outlier<-(allsum$centre < 11520.5 | allsum$centre > 20469.5)

ggplot(allsum,aes(centre,fill=outlier))+
  geom_histogram(binwidth=550)+
  theme_bw()

# Fractions:

length(which(allsum_md40$outlier == FALSE))/length(allsum_md40$outlier)
# 0.5733333
length(which(allsum_md40$centre < 11520.5))/length(allsum_md40$outlier)
# 0.26
length(which(allsum_md40$centre > 20469.5))/length(allsum_md40$outlier)
# 0.1667

## compare with snps outside centre CIs (based on HI clines)

allsum_md40<-read.csv("allsum_md40.csv")
allsum_md40$outlier<-(allsum_md40$centre < 15510 | allsum_md40$centre > 16654)
#CI:[15,510-16,654]
ggplot(allsum_md40,aes(centre,fill=outlier))+
  geom_histogram(binwidth=550)+
  theme_bw()

## Fractions:

length(which(allsum_md40$outlier == FALSE))/length(allsum_md40$outlier)
# 0.08
length(which(allsum_md40$centre < 15510))/length(allsum_md40$outlier)
# 0.4466667
length(which(allsum_md40$centre > 16654))/length(allsum_md40$outlier)
# 0.47333
