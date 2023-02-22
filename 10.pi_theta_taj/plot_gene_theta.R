library(data.table)
library(ggplot2)
library(plyr)
library(tidyverse)

#This script is used to produce theta distributions for population comparisons of Chorthippus species.
setwd("C:/Users/Jag/Documents/LMU/work/01.capture/06.pi_theta_taj/")

parcol="#A6CEED"
hybridcol="#FAD9BE"
erycol="#ED797C"
out="black"

#set colour according to rest of analyses
colperpop <- c("OUT"=out,"GAB"=parcol ,"HER"=parcol,"SOQ"=parcol,"TOU"=hybridcol,"ARA"=hybridcol,"POR"=hybridcol,"MUL"=hybridcol,"TRO"=hybridcol,"FOR"=erycol,"PAZ"=erycol,"LAN"=erycol,"ESC"=erycol,"STP"=parcol,"AIN"=parcol,"URD"=hybridcol,"OTS"=hybridcol,"ARI"=parcol,"LEK"=parcol,"MAR"=parcol,"ALM"=hybridcol,"VEN"=hybridcol)

#Read in per gene thetas estimates for all populations
getwd()
pops <- c("OUT","GAB","HER","SOQ","TOU","ARA","POR","MUL","TRO","FOR","PAZ","LAN","ESC","STP","AIN","URD","OTS","ARI","LEK","MAR","ALM","VEN")

for (pop in pops){
  sumstats<-fread(paste0("per_gene_thetas_",pop),header = TRUE)
  sumstats$nsites <- sumstats$V4 - sumstats$V3
  sumstats$pigene <- sumstats$pi/sumstats$nsites
  sumstats$thetaWgene <- sumstats$thetaW/sumstats$nsites
  assign(paste0("sumstats_",pop),sumstats)
}

thetas <- cbind.data.frame(sumstats_OUT$thetaWgene,sumstats_GAB$thetaWgene,sumstats_HER$thetaWgene,sumstats_SOQ$thetaWgene,sumstats_TOU$thetaWgene,sumstats_ARA$thetaWgene,sumstats_POR$thetaWgene,sumstats_MUL$thetaWgene,sumstats_TRO$thetaWgene,sumstats_FOR$thetaWgene,sumstats_PAZ$thetaWgene,sumstats_LAN$thetaWgene,sumstats_ESC$thetaWgene,sumstats_STP$thetaWgene,sumstats_AIN$thetaWgene,sumstats_URD$thetaWgene,sumstats_OTS$thetaWgene,sumstats_ARI$thetaWgene,sumstats_LEK$thetaWgene,sumstats_MAR$thetaWgene,sumstats_ALM$thetaWgene,sumstats_VEN$thetaWgene)
colnames(thetas)<-pops

pis <- cbind.data.frame(sumstats_OUT$pigene,sumstats_GAB$pigene,sumstats_HER$pigene,sumstats_SOQ$pigene,sumstats_TOU$pigene,sumstats_ARA$pigene,sumstats_POR$pigene,sumstats_MUL$pigene,sumstats_TRO$pigene,sumstats_FOR$pigene,sumstats_PAZ$pigene,sumstats_LAN$pigene,sumstats_ESC$pigene,sumstats_STP$pigene,sumstats_AIN$pigene,sumstats_URD$pigene,sumstats_OTS$pigene,sumstats_ARI$pigene,sumstats_LEK$pigene,sumstats_MAR$pigene,sumstats_ALM$pigene,sumstats_VEN$pigene)
colnames(pis)<-pops

tajDs <- cbind.data.frame(sumstats_OUT$TajimaD,sumstats_GAB$TajimaD,sumstats_HER$TajimaD,sumstats_SOQ$TajimaD,sumstats_TOU$TajimaD,sumstats_ARA$TajimaD,sumstats_POR$TajimaD,sumstats_MUL$TajimaD,sumstats_TRO$TajimaD,sumstats_FOR$TajimaD,sumstats_PAZ$TajimaD,sumstats_LAN$TajimaD,sumstats_ESC$TajimaD,sumstats_STP$TajimaD,sumstats_AIN$TajimaD,sumstats_URD$TajimaD,sumstats_OTS$TajimaD,sumstats_ARI$TajimaD,sumstats_LEK$TajimaD,sumstats_MAR$TajimaD,sumstats_ALM$TajimaD,sumstats_VEN$TajimaD)
colnames(tajDs)<-pops

mean(tajDs$POR)
#averagepi<- mean(c(mean(pis$POR),mean(pis$CSY),mean(pis$ERY),mean(pis$PHZ),mean(pis$PAR),mean(pis$TAR),mean(pis$AHZ),mean(pis$GOM),mean(pis$SLO),mean(pis$GOM)))
#averagepi

write.table(pis, "per_gene_pi.txt",sep = "\t", row.names = TRUE, col.names = TRUE)
write.table(thetas, "per_gene_theta.txt",sep = "\t", row.names = TRUE, col.names = TRUE)
write.table(tajDs, "per_gene_tajd.txt",sep = "\t", row.names = TRUE, col.names = TRUE)


#stack populations in same column
thetas <- stack(thetas)
colnames(thetas)<-c("theta","pop")
pis <- stack(pis)
colnames(pis) <- c("pi","pop")
tajDs <- stack(tajDs)
colnames(tajDs) <- c("tajD","pop")

#plot violin plots
pthetas <- ggplot(thetas, aes(x=theta, y=pop, fill=pop)) + 
  geom_violin(trim = TRUE, color="white" ) +
  scale_fill_manual(values=colperpop,name="Populations",
                    breaks=pops,
                    labels=c("Outgroup","Gabas","Hermine","Soques","Tourmont","Araille","Portalet","Mulas","La Troya","Formigal","Pazino","Lanuza","Escarilla","St Pee","Ainoha","Urdax","Otsondo","Arizcun","Lekaroz","Mardea","Almandoz","Ventas de Arraiz")) +
  stat_summary(fun=mean, geom="point", size=2, color = "white") +
  #stat_summary(fun=median, geom="point", shape=21, size=2, color="white")+
  scale_y_discrete(limits=rev) +
  xlim(0,0.0225)
pthetas

ppis <- ggplot(pis, aes(x=pi, y=pop, fill=pop)) + 
  geom_violin(trim = TRUE, color="white" ) +
  scale_fill_manual(values=colperpop,name="Populations",
                    breaks=pops,
                    labels=c("Outgroup","Gabas","Hermine","Soques","Tourmont","Araille","Portalet","Mulas","La Troya","Formigal","Pazino","Lanuza","Escarilla","St Pee","Ainoha","Urdax","Otsondo","Arizcun","Lekaroz","Mardea","Almandoz","Ventas de Arraiz")) +
  stat_summary(fun=mean, geom="point", size=2, color = "white") +
  #stat_summary(fun=median, geom="point", shape=21, size=2, color="white")+
  scale_y_discrete(limits=rev) +
  xlim(0,0.0225)

ppis


ptajDs <- ggplot(tajDs, aes(x=tajD, y=pop, fill=pop)) + 
  geom_violin(trim = TRUE, color="white" ) +
  scale_fill_manual(values=colperpop,name="Populations",
                    breaks=pops,
                    labels=c("Outgroup","Gabas","Hermine","Soques","Tourmont","Araille","Portalet","Mulas","La Troya","Formigal","Pazino","Lanuza","Escarilla","St Pee","Ainoha","Urdax","Otsondo","Arizcun","Lekaroz","Mardea","Almandoz","Ventas de Arraiz")) +
  stat_summary(fun=mean, geom="point", size=2, color = "white") +
  #stat_summary(fun=median, geom="point", shape=21, size=2, color="white")+
  scale_y_discrete(limits=rev) + 
  xlim(-2,2)
  #coord_fixed(ratio = 0.5, xlim = NULL, ylim = NULL, expand = TRUE, clip = "on")

ptajDs

ggsave("sizedtajD",ptajDs)

