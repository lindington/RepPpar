#####################################################
#####################################################
######### SNP SELECTION FOR CAPTURE CLINES  ######### 
#####################################################
#####################################################

###########################################
######### LOAD PROGRAMS AND FILES ######### 
###########################################

### load packages
library("dplyr")
library("tidyverse")
library("ggplot2")
library("data.table")

### set to working directory
getwd()
setwd("C:/Users/Jag/Documents/LMU/work/01.capture/07.hzar")

# name cols of mafs files
colns <-c("chromo","position","major","minor","ref","anc","knownEM","pK-EM","nInd")

# average allele freq of extremes: PAR includes GABAS & HERMINE, ERY includes LANUZA & ESCARILLA
pops <- c("PAR", "ERY")

## load linked allele freq files
PAR_linked <- read.csv("af_cline_PAR.mafs", sep="\t", header=TRUE, row.names = NULL, col.names = colns, colClasses = c(rep(NA,5),"NULL",NA,"NULL","NULL"))
ERY_linked <- read.csv("af_cline_ERY.mafs", sep="\t", header=TRUE, row.names = NULL, col.names = colns, colClasses = c(rep(NA,5),"NULL",NA,"NULL","NULL"))

## load locus information file "bait.info"
loci<-read.csv("bait.info",sep="\t",header=FALSE,row.names = NULL,colClasses = c(rep(NA,4),rep("NULL",4)),col.names = c("CHR","BEG","END","LOCUS",rep("x",4)))

## add locus inf to allele freq files by replacing redundant "chromo" with locus info)
for (i in 1:length(loci$CHR)){
  ks <- (which(ERY_linked$position %in% loci[[i,2]]:loci[[i,3]]));
  for (k in ks){
    ERY_linked$chromo[k]<-loci$LOCUS[i]
  }
}

for (i in 1:length(loci$CHR)){
  ks <- (which(PAR_linked$position %in% loci[[i,2]]:loci[[i,3]]));
  for (k in ks){
    PAR_linked$chromo[k]<-loci$LOCUS[i]
  }
}

## remove loop placeholders
remove(i,ks,k)

######################################################################################
######## Calculate difference between allele frequencies in transect extremes ########
######################################################################################

colnames(ais40md_l)<-c("CHROM","POS","MAJ_ALLELE","MAJ_FREQ","MIN_ALLELE","MIN_FREQ")
colnames(bar40md_l)<-c("CHROM","POS","MAJ_ALLELE","MAJ_FREQ","MIN_ALLELE","MIN_FREQ")

ais40md_l <- ais40md_l %>%
  unite('CHROM_POS',CHROM:POS, remove = TRUE)
bar40md_l <- bar40md_l %>%
  unite('CHROM_POS',CHROM:POS, remove = TRUE)

### make dataframe with positions found in both extremes
both_linked <- as.data.frame(inner_join(ERY_linked, PAR_linked, by="position"))
both_linked <- data.frame(both_linked$chromo.x, both_linked$position, both_linked$ref.x, both_linked$major.x, as.numeric(both_linked$knownEM.x), both_linked$major.y, as.numeric(both_linked$knownEM.y), both_linked$minor.x, as.numeric(0.5-both_linked$knownEM.x), both_linked$minor.y, as.numeric(0.5-both_linked$knownEM.y))

colnames(both_linked)<-c("CHROM","POS","REF","MAJ_ALLELE_ERY","MAJ_FREQ_ERY","MAJ_ALLELE_PAR","MAJ_FREQ_PAR","MIN_ALLELE_ERY","MIN_FREQ_ERY","MIN_ALLELE_PAR","MIN_FREQ_PAR")

haf_alp40md_l <- data.frame(alp40md_l$CHROM_POS,alp40md_l$MAJ_ALLELE_AIS,alp40md_l$MAJ_FREQ_AIS-alp40md_l$MAJ_FREQ_BAR, abs(alp40md_l$MAJ_FREQ_AIS-alp40md_l$MAJ_FREQ_BAR))
colnames(haf_alp40md_l) <- c("CHROM_POS","ALLELE","DELTA_FREQ","DELTA_FREQ_ABS")

##making sure the major allele is the same (should be, because they're from the same assembly)
all(ERY_linked$major==PAR_linked$major)

##removing SNPs that have a frequency change of 0
haf_alp40md_l_onlypos<-haf_alp40md_l[haf_alp40md_l$DELTA_FREQ != 0, ]

###plot distribution of delta af

ggplot(haf_alp40md_l_onlypos, aes(DELTA_FREQ_ABS)) +
  geom_density(alpha = 0.1,adjust = 0.01, colour="darkred", fill="red")

#############################################################
##################### SAME FOR PYRENEES ##################### 
#############################################################


gab40md_l <- read.csv("af_GAB_md40_l.frq", sep="\t", header=TRUE, row.names = NULL, col.names = colns)
esc40md_l <- read.csv("af_ESC_md40_l.frq", sep="\t", header=TRUE, row.names = NULL, col.names = colns)
gab40md_l.info1 <- data.frame(t(sapply(as.character(gab40md_l[,5]), function(y) strsplit(y,split=":")[[1]])))
gab40md_l.info2 <- data.frame(t(sapply(as.character(gab40md_l[,6]), function(y) strsplit(y,split=":")[[1]])))
esc40md_l.info1 <- data.frame(t(sapply(as.character(esc40md_l[,5]), function(y) strsplit(y,split=":")[[1]])))
esc40md_l.info2 <- data.frame(t(sapply(as.character(esc40md_l[,6]), function(y) strsplit(y,split=":")[[1]])))
colnames(gab40md_l.info1)<-c("MAJ_ALLELE","MAJ_FREQ")
colnames(gab40md_l.info2)<-c("MIN_ALLELE","MIN_FREQ")
colnames(esc40md_l.info1)<-c("MAJ_ALLELE","MAJ_FREQ")
colnames(esc40md_l.info2)<-c("MIN_ALLELE","MIN_FREQ")

## dont include N_ALLELES because always =2
#gab40md_l$N_ALLELES[!(gab40md_l$N_ALLELES %in% c('2'))]

gab40md_l <- data.frame(gab40md_l$CHROM, gab40md_l$POS, gab40md_l$N_CHR, gab40md_l.info1$MAJ_ALLELE, gab40md_l.info1$MAJ_FREQ, gab40md_l.info2$MIN_ALLELE, gab40md_l.info2$MIN_FREQ)
esc40md_l <- data.frame(esc40md_l$CHROM, esc40md_l$POS, esc40md_l$N_CHR, esc40md_l.info1$MAJ_ALLELE, esc40md_l.info1$MAJ_FREQ, esc40md_l.info2$MIN_ALLELE, esc40md_l.info2$MIN_FREQ)
remove(gab40md_l.info1,gab40md_l.info2,esc40md_l.info1,esc40md_l.info2)

colnames(gab40md_l)<-c("CHROM","POS","N_CHR","MAJ_ALLELE","MAJ_FREQ","MIN_ALLELE","MIN_FREQ")
colnames(esc40md_l)<-c("CHROM","POS","N_CHR","MAJ_ALLELE","MAJ_FREQ","MIN_ALLELE","MIN_FREQ")

gab40md_l <- gab40md_l %>%
  unite('CHROM_POS',CHROM:POS, remove = TRUE)
esc40md_l <- esc40md_l %>%
  unite('CHROM_POS',CHROM:POS, remove = TRUE)

### make dataframe with positions found in both extremes
pyr40md_l <- as.data.frame(inner_join(esc40md_l, gab40md_l, by="CHROM_POS"))
colnames(pyr40md_l)<-c("CHROM_POS","N_CHR_gab","MAJ_ALLELE_gab","MAJ_FREQ_gab","MIN_ALLELE_gab","MIN_FREQ_gab","N_CHR_esc","MAJ_ALLELE_esc","MAJ_FREQ_esc","MIN_ALLELE_esc","MIN_FREQ_esc")
pyr40md_l$MAJ_FREQ_gab<-as.numeric(pyr40md_l$MAJ_FREQ_gab)
pyr40md_l$MIN_FREQ_gab<-as.numeric(pyr40md_l$MIN_FREQ_gab)
pyr40md_l$MAJ_FREQ_esc<-as.numeric(pyr40md_l$MAJ_FREQ_esc)
pyr40md_l$MIN_FREQ_esc<-as.numeric(pyr40md_l$MIN_FREQ_esc)


haf_pyr40md_l <- data.frame(pyr40md_l$CHROM_POS,pyr40md_l$MAJ_ALLELE_gab,pyr40md_l$MAJ_FREQ_gab-pyr40md_l$MAJ_FREQ_esc, abs(pyr40md_l$MAJ_FREQ_gab-pyr40md_l$MAJ_FREQ_esc))
colnames(haf_pyr40md_l) <- c("CHROM_POS","ALLELE","DELTA_FREQ","DELTA_FREQ_ABS")

##making sure the major allele is the same
all(pyr40md_l$MAJ_ALLELE_gab==haf_pyr40md_l$ALLELE)

haf_pyr40md_l_onlypos<-haf_pyr40md_l[haf_pyr40md_l$DELTA_FREQ != 0, ]

###plot distribution of delta af

ggplot(haf_pyr40md_l_onlypos, aes(DELTA_FREQ_ABS)) +
  geom_density(alpha = 0.1,adjust = 0.01, colour="darkred", fill="red")
#+  xlim(-0.25, 0.25)



#############################################################
#############################################################
############# Get highest delta SNP per locus ############### 
#############################################################
#############################################################



popsalp <- c("AIS","PIA", "SAM","PON", "COL", "CER","CON", "BAR")
popspyr <- c("GAB","HER", "SOQ", "TOU", "ARA", "POR", "MUL", "FOR", "PAZ","LAN", "ESC")

colns <-c("CHROM","POS","N_ALLELES","N_CHR","MAJ_ALLELE_FREQ","MIN_ALLELE_FREQ") 

alp.list <- list(read.csv("af_AIS_md40_l.frq", sep="\t", header=TRUE, row.names =  NULL, col.names = colns),
                 read.csv("af_PIA_md40_l.frq", sep="\t", header=TRUE, row.names =  NULL, col.names = colns),
                 read.csv("af_SAM_md40_l.frq", sep="\t", header=TRUE, row.names =  NULL, col.names = colns), 
                 read.csv("af_PON_md40_l.frq", sep="\t", header=TRUE, row.names =  NULL, col.names = colns), 
                 read.csv("af_COL_md40_l.frq", sep="\t", header=TRUE, row.names =  NULL, col.names = colns), 
                 read.csv("af_CER_md40_l.frq", sep="\t", header=TRUE, row.names =  NULL, col.names = colns), 
                 read.csv("af_CON_md40_l.frq", sep="\t", header=TRUE, row.names =  NULL, col.names = colns),
                 read.csv("af_BAR_md40_l.frq", sep="\t", header=TRUE, row.names =  NULL, col.names = colns))


pyr.list <- list(read.csv("af_GAB_md40_l.frq", sep="\t", header=TRUE, row.names =  NULL, col.names = colns), 
                 read.csv("af_HER_md40_l.frq", sep="\t", header=TRUE, row.names =  NULL, col.names = colns), 
                 read.csv("af_SOQ_md40_l.frq", sep="\t", header=TRUE, row.names =  NULL, col.names = colns), 
                 read.csv("af_TOU_md40_l.frq", sep="\t", header=TRUE, row.names =  NULL, col.names = colns), 
                 read.csv("af_ARA_md40_l.frq", sep="\t", header=TRUE, row.names =  NULL, col.names = colns), 
                 read.csv("af_POR_md40_l.frq", sep="\t", header=TRUE, row.names =  NULL, col.names = colns), 
                 read.csv("af_MUL_md40_l.frq", sep="\t", header=TRUE, row.names =  NULL, col.names = colns), 
                 read.csv("af_FOR_md40_l.frq", sep="\t", header=TRUE, row.names =  NULL, col.names = colns), 
                 read.csv("af_PAZ_md40_l.frq", sep="\t", header=TRUE, row.names =  NULL, col.names = colns), 
                 read.csv("af_LAN_md40_l.frq", sep="\t", header=TRUE, row.names =  NULL, col.names = colns), 
                 read.csv("af_ESC_md40_l.frq", sep="\t", header=TRUE, row.names =  NULL, col.names = colns))

## lapply applies functions to every element of a list. 
alp.list <- lapply(alp.list, function(x) x %>% 
                     unite('CHROM_POS',CHROM:POS, remove = TRUE) %>% 
                     separate(MAJ_ALLELE_FREQ, c("MAJ_ALLELE","MAJ_FREQ"), sep=":", remove = TRUE) %>% 
                     separate(MIN_ALLELE_FREQ, c("MIN_ALLELE","MIN_FREQ"), sep=":", remove = TRUE) %>% 
                     filter(CHROM_POS %in% haf_alp40md_l_onlypos$CHROM_POS))

length(alp.list[[4]]$CHROM_POS)

snps_alp40md <- haf_alp40md_l_onlypos %>% 
  filter(CHROM_POS %in% alp.list[[1]]$CHROM_POS) %>%
  filter(CHROM_POS %in% alp.list[[2]]$CHROM_POS) %>%
  filter(CHROM_POS %in% alp.list[[3]]$CHROM_POS) %>%
  filter(CHROM_POS %in% alp.list[[4]]$CHROM_POS) %>%
  filter(CHROM_POS %in% alp.list[[5]]$CHROM_POS) %>%
  filter(CHROM_POS %in% alp.list[[6]]$CHROM_POS) %>%
  filter(CHROM_POS %in% alp.list[[7]]$CHROM_POS) %>%
  filter(CHROM_POS %in% alp.list[[8]]$CHROM_POS)

#making sure all snps are in all pops
#all(snps_alp40md$CHROM_POS %in% alp.list[[6]]$CHROM_POS)

##find highest delta freq abs
hl_alp40md <- separate(snps_alp40md, CHROM_POS, c("SEQ","CHROM", "POS"), sep="_")
hl_alp40md <-hl_alp40md %>% unite('CHROM',SEQ:CHROM, remove = TRUE)
hl_alp40md <- hl_alp40md %>% group_by(CHROM) %>% slice(which.max(DELTA_FREQ_ABS))
hl_alp40md <- hl_alp40md %>% unite('CHROM_POS',CHROM:POS, remove = TRUE)

#final list of cline snps:
clineSNPs_alp40md <- hl_alp40md$CHROM_POS
fwrite(list(clineSNPs_alp40md), file="../clineSNPs_alp40md.txt")

##filter all pops by highest snp (present in all pops) per locus.
alp.list <- lapply(alp.list, function(x) x %>% filter(CHROM_POS %in% clineSNPs_alp40md))


alp_frqs_md40 <- data.frame(alp.list[[1]]$MAJ_FREQ, alp.list[[2]]$MAJ_FREQ, alp.list[[3]]$MAJ_FREQ, alp.list[[4]]$MAJ_FREQ, alp.list[[5]]$MAJ_FREQ, alp.list[[6]]$MAJ_FREQ, alp.list[[7]]$MAJ_FREQ, alp.list[[8]]$MAJ_FREQ)
alp_frqs_md40 <- as.data.frame(t(alp_frqs_md40))
colnames(alp_frqs_md40)<-alp.list[[1]]$CHROM_POS
rownames(alp_frqs_md40)<-NULL

alpsamsis_md40 <- data.frame(alp.list[[1]]$N_CHR, alp.list[[2]]$N_CHR, alp.list[[3]]$N_CHR, alp.list[[4]]$N_CHR, alp.list[[5]]$N_CHR, alp.list[[6]]$N_CHR, alp.list[[7]]$N_CHR, alp.list[[8]]$N_CHR)
alpsamsis_md40 <- as.data.frame(t(alpsamsis_md40))
colnames(alpsamsis_md40)<-paste("N_",alp.list[[1]]$CHROM_POS,sep = "")
rownames(alpsamsis_md40)<-NULL

##alternate allelefrequency and allele sample size columns

alp_frqs_md40 <-cbind.data.frame(alp_frqs_md40,alpsamsis_md40)
alp_frqs_md40 <- as.data.frame(alp_frqs_md40[, c(matrix(1:ncol(alp_frqs_md40), nrow = 2, byrow = T))])

alp_frqs_md40$distance<-alpdist
alp_frqs_md40$population<-popsalp
alp_frqs_md40 <- alp_frqs_md40 %>%
  relocate(distance) %>%
  relocate(population)
write.csv(alp_frqs_md40,"../alp_frqs_md40.csv", row.names = FALSE)

## lapply applies functions to every element of a list. 
pyr.list <- lapply(pyr.list, function(x) x %>% 
                     unite('CHROM_POS',CHROM:POS, remove = TRUE) %>% 
                     separate(MAJ_ALLELE_FREQ, c("MAJ_ALLELE","MAJ_FREQ"), sep=":", remove = TRUE) %>% 
                     separate(MIN_ALLELE_FREQ, c("MIN_ALLELE","MIN_FREQ"), sep=":", remove = TRUE) %>% 
                     filter(CHROM_POS %in% haf_pyr40md_l_onlypos$CHROM_POS))

snps_pyr40md <- haf_pyr40md_l_onlypos %>% 
  filter(CHROM_POS %in% pyr.list[[1]]$CHROM_POS) %>%
  filter(CHROM_POS %in% pyr.list[[2]]$CHROM_POS) %>%
  filter(CHROM_POS %in% pyr.list[[3]]$CHROM_POS) %>%
  filter(CHROM_POS %in% pyr.list[[4]]$CHROM_POS) %>%
  filter(CHROM_POS %in% pyr.list[[5]]$CHROM_POS) %>%
  filter(CHROM_POS %in% pyr.list[[6]]$CHROM_POS) %>%
  filter(CHROM_POS %in% pyr.list[[7]]$CHROM_POS) %>%
  filter(CHROM_POS %in% pyr.list[[8]]$CHROM_POS) %>%
  filter(CHROM_POS %in% pyr.list[[9]]$CHROM_POS) %>%
  filter(CHROM_POS %in% pyr.list[[10]]$CHROM_POS) %>%
  filter(CHROM_POS %in% pyr.list[[11]]$CHROM_POS)
  
#making sure all snps are in all pops
#all(snps_pyr40md$CHROM_POS %in% pyr.list[[6]]$CHROM_POS)

##find highest delta freq abs
hl_pyr40md <- separate(snps_pyr40md, CHROM_POS, c("SEQ","CHROM", "POS"), sep="_")
hl_pyr40md <-hl_pyr40md %>% unite('CHROM',SEQ:CHROM, remove = TRUE)
hl_pyr40md <- hl_pyr40md %>% group_by(CHROM) %>% slice(which.max(DELTA_FREQ_ABS))
hl_pyr40md <- hl_pyr40md %>% unite('CHROM_POS',CHROM:POS, remove = TRUE)

#final list of cline snps:
clineSNPs_pyr40md <- hl_pyr40md$CHROM_POS
fwrite(list(clineSNPs_pyr40md), file="../clineSNPs_pyr40md.txt")

##filter all pops by highest snp (present in all pops) per locus.
pyr.list <- lapply(pyr.list, function(x) x %>% filter(CHROM_POS %in% clineSNPs_pyr40md))

pyr_frqs_md40 <- data.frame(pyr.list[[1]]$MAJ_FREQ, pyr.list[[2]]$MAJ_FREQ, pyr.list[[3]]$MAJ_FREQ, pyr.list[[4]]$MAJ_FREQ, pyr.list[[5]]$MAJ_FREQ, pyr.list[[6]]$MAJ_FREQ, pyr.list[[7]]$MAJ_FREQ, pyr.list[[8]]$MAJ_FREQ, pyr.list[[9]]$MAJ_FREQ, pyr.list[[10]]$MAJ_FREQ, pyr.list[[11]]$MAJ_FREQ)
pyr_frqs_md40 <- as.data.frame(t(pyr_frqs_md40))
colnames(pyr_frqs_md40)<-pyr.list[[1]]$CHROM_POS
rownames(pyr_frqs_md40)<-NULL

pyrsamsis_md40 <- data.frame(pyr.list[[1]]$N_CHR, pyr.list[[2]]$N_CHR, pyr.list[[3]]$N_CHR, pyr.list[[4]]$N_CHR, pyr.list[[5]]$N_CHR, pyr.list[[6]]$N_CHR, pyr.list[[7]]$N_CHR, pyr.list[[8]]$N_CHR, pyr.list[[9]]$N_CHR, pyr.list[[10]]$N_CHR, pyr.list[[11]]$N_CHR)
pyrsamsis_md40 <- as.data.frame(t(pyrsamsis_md40))
colnames(pyrsamsis_md40)<-paste("N_",pyr.list[[1]]$CHROM_POS,sep = "")
rownames(pyrsamsis_md40)<-NULL

##alternate allelefrequency and allele sample size columns

pyr_frqs_md40 <-cbind.data.frame(pyr_frqs_md40,pyrsamsis_md40)
pyr_frqs_md40 <- as.data.frame(pyr_frqs_md40[, c(matrix(1:ncol(pyr_frqs_md40), nrow = 2, byrow = T))])

pyr_frqs_md40$distance<-pyrdist
pyr_frqs_md40$population<-popspyr
pyr_frqs_md40 <- pyr_frqs_md40 %>%
  relocate(distance) %>%
  relocate(population)
write.csv(pyr_frqs_md40,"../pyr_frqs_md40.csv", row.names = FALSE)


##orientation of clines may be relevant -> invert upwards

#curious about overlap between highest snp, but prob overfiltered
#md40snp<- intersect(clineSNPs_alp40md, clineSNPs_pyr40md)
#md60snp<- intersect(clineSNPs_alp60md, clineSNPs_pyr60md)
#mdbothsnp<- intersect(md40snp, md60snp)
######################################################################################
######################################################################################
####################################### 60 md ########################################
######################################################################################
######################################################################################

######################################################################################
######## Calculate difference between allele frequencies in transect extremes ########
######################################################################################

ais60md_l <- read.csv("af_ais_md60_l.frq", sep="\t", header=TRUE, row.names = NULL, col.names = colns)
bar60md_l <- read.csv("af_bar_md60_l.frq", sep="\t", header=TRUE, row.names = NULL, col.names = colns)
ais60md_l.info1 <- data.frame(t(sapply(as.character(ais60md_l[,5]), function(y) strsplit(y,split=":")[[1]])))
ais60md_l.info2 <- data.frame(t(sapply(as.character(ais60md_l[,6]), function(y) strsplit(y,split=":")[[1]])))
bar60md_l.info1 <- data.frame(t(sapply(as.character(bar60md_l[,5]), function(y) strsplit(y,split=":")[[1]])))
bar60md_l.info2 <- data.frame(t(sapply(as.character(bar60md_l[,6]), function(y) strsplit(y,split=":")[[1]])))
colnames(ais60md_l.info1)<-c("MAJ_ALLELE","MAJ_FREQ")
colnames(ais60md_l.info2)<-c("MIN_ALLELE","MIN_FREQ")
colnames(bar60md_l.info1)<-c("MAJ_ALLELE","MAJ_FREQ")
colnames(bar60md_l.info2)<-c("MIN_ALLELE","MIN_FREQ")

## dont include N_ALLELES because always =2
#ais60md_l$N_ALLELES[!(ais60md_l$N_ALLELES %in% c('2'))]

ais60md_l <- data.frame(ais60md_l$CHROM, ais60md_l$POS, ais60md_l$N_CHR, ais60md_l.info1$MAJ_ALLELE, ais60md_l.info1$MAJ_FREQ, ais60md_l.info2$MIN_ALLELE, ais60md_l.info2$MIN_FREQ)
bar60md_l <- data.frame(bar60md_l$CHROM, bar60md_l$POS, bar60md_l$N_CHR, bar60md_l.info1$MAJ_ALLELE, bar60md_l.info1$MAJ_FREQ, bar60md_l.info2$MIN_ALLELE, bar60md_l.info2$MIN_FREQ)
remove(ais60md_l.info1,ais60md_l.info2,bar60md_l.info1,bar60md_l.info2)

colnames(ais60md_l)<-c("CHROM","POS","N_CHR","MAJ_ALLELE","MAJ_FREQ","MIN_ALLELE","MIN_FREQ")
colnames(bar60md_l)<-c("CHROM","POS","N_CHR","MAJ_ALLELE","MAJ_FREQ","MIN_ALLELE","MIN_FREQ")

ais60md_l <- ais60md_l %>%
  unite('CHROM_POS',CHROM:POS, remove = TRUE)
bar60md_l <- bar60md_l %>%
  unite('CHROM_POS',CHROM:POS, remove = TRUE)

### make dataframe with positions found in both extremes
alp60md_l <- as.data.frame(inner_join(bar60md_l, ais60md_l, by="CHROM_POS"))
colnames(alp60md_l)<-c("CHROM_POS","N_CHR_AIS","MAJ_ALLELE_AIS","MAJ_FREQ_AIS","MIN_ALLELE_AIS","MIN_FREQ_AIS","N_CHR_BAR","MAJ_ALLELE_BAR","MAJ_FREQ_BAR","MIN_ALLELE_BAR","MIN_FREQ_BAR")
alp60md_l$MAJ_FREQ_AIS<-as.numeric(alp60md_l$MAJ_FREQ_AIS)
alp60md_l$MIN_FREQ_AIS<-as.numeric(alp60md_l$MIN_FREQ_AIS)
alp60md_l$MAJ_FREQ_BAR<-as.numeric(alp60md_l$MAJ_FREQ_BAR)
alp60md_l$MIN_FREQ_BAR<-as.numeric(alp60md_l$MIN_FREQ_BAR)

haf_alp60md_l <- data.frame(alp60md_l$CHROM_POS,alp60md_l$MAJ_ALLELE_AIS,alp60md_l$MAJ_FREQ_AIS-alp60md_l$MAJ_FREQ_BAR, abs(alp60md_l$MAJ_FREQ_AIS-alp60md_l$MAJ_FREQ_BAR))
colnames(haf_alp60md_l) <- c("CHROM_POS","ALLELE","DELTA_FREQ","DELTA_FREQ_ABS")

##making sure the major allele is the same (should be, because they're from the same assembly)
all(alp60md_l$MAJ_ALLELE_AIS==haf_alp60md_l$ALLELE)

##removing SNPs that have a frequency change of 0
haf_alp60md_l_onlypos<-haf_alp60md_l[haf_alp60md_l$DELTA_FREQ != 0, ]

###plot distribution of delta af

hr <- ggplot(haf_alp60md_l_onlypos, aes(DELTA_FREQ_ABS)) + 
  geom_density(alpha = 0.1,adjust = 0.01, colour="darkred", fill="red")
lr <- ggplot(haf_alp60md_l_onlypos, aes(DELTA_FREQ_ABS)) + 
  geom_density(alpha = 0.1,adjust = 0.5, colour="darkred", fill="red")

#############################################################
##################### SAME FOR PYRENEES ##################### 
#############################################################


gab60md_l <- read.csv("af_GAB_md60_l.frq", sep="\t", header=TRUE, row.names = NULL, col.names = colns)
esc60md_l <- read.csv("af_ESC_md60_l.frq", sep="\t", header=TRUE, row.names = NULL, col.names = colns)
gab60md_l.info1 <- data.frame(t(sapply(as.character(gab60md_l[,5]), function(y) strsplit(y,split=":")[[1]])))
gab60md_l.info2 <- data.frame(t(sapply(as.character(gab60md_l[,6]), function(y) strsplit(y,split=":")[[1]])))
esc60md_l.info1 <- data.frame(t(sapply(as.character(esc60md_l[,5]), function(y) strsplit(y,split=":")[[1]])))
esc60md_l.info2 <- data.frame(t(sapply(as.character(esc60md_l[,6]), function(y) strsplit(y,split=":")[[1]])))
colnames(gab60md_l.info1)<-c("MAJ_ALLELE","MAJ_FREQ")
colnames(gab60md_l.info2)<-c("MIN_ALLELE","MIN_FREQ")
colnames(esc60md_l.info1)<-c("MAJ_ALLELE","MAJ_FREQ")
colnames(esc60md_l.info2)<-c("MIN_ALLELE","MIN_FREQ")

## dont include N_ALLELES because always =2
#gab60md_l$N_ALLELES[!(gab60md_l$N_ALLELES %in% c('2'))]

gab60md_l <- data.frame(gab60md_l$CHROM, gab60md_l$POS, gab60md_l$N_CHR, gab60md_l.info1$MAJ_ALLELE, gab60md_l.info1$MAJ_FREQ, gab60md_l.info2$MIN_ALLELE, gab60md_l.info2$MIN_FREQ)
esc60md_l <- data.frame(esc60md_l$CHROM, esc60md_l$POS, esc60md_l$N_CHR, esc60md_l.info1$MAJ_ALLELE, esc60md_l.info1$MAJ_FREQ, esc60md_l.info2$MIN_ALLELE, esc60md_l.info2$MIN_FREQ)
remove(gab60md_l.info1,gab60md_l.info2,esc60md_l.info1,esc60md_l.info2)

colnames(gab60md_l)<-c("CHROM","POS","N_CHR","MAJ_ALLELE","MAJ_FREQ","MIN_ALLELE","MIN_FREQ")
colnames(esc60md_l)<-c("CHROM","POS","N_CHR","MAJ_ALLELE","MAJ_FREQ","MIN_ALLELE","MIN_FREQ")

gab60md_l <- gab60md_l %>%
  unite('CHROM_POS',CHROM:POS, remove = TRUE)
esc60md_l <- esc60md_l %>%
  unite('CHROM_POS',CHROM:POS, remove = TRUE)

### make dataframe with positions found in both extremes
pyr60md_l <- as.data.frame(inner_join(esc60md_l, gab60md_l, by="CHROM_POS"))
colnames(pyr60md_l)<-c("CHROM_POS","N_CHR_gab","MAJ_ALLELE_gab","MAJ_FREQ_gab","MIN_ALLELE_gab","MIN_FREQ_gab","N_CHR_esc","MAJ_ALLELE_esc","MAJ_FREQ_esc","MIN_ALLELE_esc","MIN_FREQ_esc")
pyr60md_l$MAJ_FREQ_gab<-as.numeric(pyr60md_l$MAJ_FREQ_gab)
pyr60md_l$MIN_FREQ_gab<-as.numeric(pyr60md_l$MIN_FREQ_gab)
pyr60md_l$MAJ_FREQ_esc<-as.numeric(pyr60md_l$MAJ_FREQ_esc)
pyr60md_l$MIN_FREQ_esc<-as.numeric(pyr60md_l$MIN_FREQ_esc)

haf_pyr60md_l <- data.frame(pyr60md_l$CHROM_POS,pyr60md_l$MAJ_ALLELE_gab,pyr60md_l$MAJ_FREQ_gab-pyr60md_l$MAJ_FREQ_esc, abs(pyr60md_l$MAJ_FREQ_gab-pyr60md_l$MAJ_FREQ_esc))
colnames(haf_pyr60md_l) <- c("CHROM_POS","ALLELE","DELTA_FREQ","DELTA_FREQ_ABS")

##making sure the major allele is the same
all(pyr60md_l$MAJ_ALLELE_gab==haf_pyr60md_l$ALLELE)

##removing SNPs that have a frequency change of 0
haf_pyr60md_l_onlypos<-haf_pyr60md_l[haf_pyr60md_l$DELTA_FREQ != 0, ]

###plot distribution of delta af

hr <- ggplot(haf_pyr60md_l_onlypos, aes(DELTA_FREQ_ABS)) +
  geom_density(alpha = 0.1,adjust = 0.01, colour="darkred", fill="red")
lr <- ggplot(haf_pyr60md_l_onlypos, aes(DELTA_FREQ_ABS)) +
  geom_density(alpha = 0.1,adjust = 0.5, colour="darkred", fill="red")

#############################################################
############# Get highest delta SNP per locus ############### 
#############################################################

## separate the locus from locus position to group by locus
### ALPS
hl_alp60md <- separate(haf_alp60md_l_onlypos, CHROM_POS, c("SEQ","CHROM", "POS"), sep="_")
hl_alp60md <-hl_alp60md %>% unite('CHROM',SEQ:CHROM, remove = TRUE)

##find highest delta freq abs
hl_alp60md <- hl_alp60md %>% group_by(CHROM) %>% slice(which.max(DELTA_FREQ_ABS))
hl_alp60md <- hl_alp60md %>% unite('CHROM_POS',CHROM:POS, remove = TRUE)

clineSNPs_alp60md <- hl_alp60md$CHROM_POS
fwrite(list(clineSNPs_alp60md), file="../clineSNPs_alp60md.txt")

### PYRS
hl_pyr60md <- separate(haf_pyr60md_l_onlypos, CHROM_POS, c("SEQ","CHROM", "POS"), sep="_")
hl_pyr60md <-hl_pyr60md %>% unite('CHROM',SEQ:CHROM, remove = TRUE)

##find highest delta freq abs 
hl_pyr60md <- hl_pyr60md %>% group_by(CHROM) %>% slice(which.max(DELTA_FREQ_ABS))
hl_pyr60md <- hl_pyr60md %>% unite('CHROM_POS',CHROM:POS, remove = TRUE)

clineSNPs_pyr60md <- hl_pyr60md$CHROM_POS
fwrite(list(clineSNPs_pyr60md), file="../clineSNPs_pyr60md.txt")

#############################################################
#############################################################
######## All pops: Get highest delta SNP per locus ########## 
#############################################################
#############################################################


alp.list <- list(read.csv("af_AIS_md60_l.frq", sep="\t", header=TRUE, row.names =  NULL, col.names = colns),
                 read.csv("af_PIA_md60_l.frq", sep="\t", header=TRUE, row.names =  NULL, col.names = colns), 
                 read.csv("af_SAM_md60_l.frq", sep="\t", header=TRUE, row.names =  NULL, col.names = colns), 
                 read.csv("af_PON_md60_l.frq", sep="\t", header=TRUE, row.names =  NULL, col.names = colns), 
                 read.csv("af_COL_md60_l.frq", sep="\t", header=TRUE, row.names =  NULL, col.names = colns), 
                 read.csv("af_CER_md60_l.frq", sep="\t", header=TRUE, row.names =  NULL, col.names = colns),
                 read.csv("af_CON_md60_l.frq", sep="\t", header=TRUE, row.names =  NULL, col.names = colns), 
                 read.csv("af_BAR_md60_l.frq", sep="\t", header=TRUE, row.names =  NULL, col.names = colns))
                 

pyr.list <- list(read.csv("af_GAB_md60_l.frq", sep="\t", header=TRUE, row.names =  NULL, col.names = colns), 
                 read.csv("af_HER_md60_l.frq", sep="\t", header=TRUE, row.names =  NULL, col.names = colns), 
                 read.csv("af_SOQ_md60_l.frq", sep="\t", header=TRUE, row.names =  NULL, col.names = colns), 
                 read.csv("af_TOU_md60_l.frq", sep="\t", header=TRUE, row.names =  NULL, col.names = colns), 
                 read.csv("af_ARA_md60_l.frq", sep="\t", header=TRUE, row.names =  NULL, col.names = colns), 
                 read.csv("af_POR_md60_l.frq", sep="\t", header=TRUE, row.names =  NULL, col.names = colns), 
                 read.csv("af_MUL_md60_l.frq", sep="\t", header=TRUE, row.names =  NULL, col.names = colns), 
                 read.csv("af_FOR_md60_l.frq", sep="\t", header=TRUE, row.names =  NULL, col.names = colns), 
                 read.csv("af_PAZ_md60_l.frq", sep="\t", header=TRUE, row.names =  NULL, col.names = colns), 
                 read.csv("af_LAN_md60_l.frq", sep="\t", header=TRUE, row.names =  NULL, col.names = colns), 
                 read.csv("af_ESC_md60_l.frq", sep="\t", header=TRUE, row.names =  NULL, col.names = colns))
   
## lapply applies functions to every element of a list. 
alp.list <- lapply(alp.list, function(x) x %>% 
                     unite('CHROM_POS',CHROM:POS, remove = TRUE) %>% 
                     separate(MAJ_ALLELE_FREQ, c("MAJ_ALLELE","MAJ_FREQ"), sep=":", remove = TRUE) %>% 
                     separate(MIN_ALLELE_FREQ, c("MIN_ALLELE","MIN_FREQ"), sep=":", remove = TRUE) %>% 
                     filter(CHROM_POS %in% haf_alp60md_l_onlypos$CHROM_POS))

length(alp.list[[4]]$CHROM_POS)

snps_alp60md <- haf_alp60md_l_onlypos %>% 
  filter(CHROM_POS %in% alp.list[[1]]$CHROM_POS) %>%
  filter(CHROM_POS %in% alp.list[[2]]$CHROM_POS) %>%
  filter(CHROM_POS %in% alp.list[[3]]$CHROM_POS) %>%
  filter(CHROM_POS %in% alp.list[[4]]$CHROM_POS) %>%
  filter(CHROM_POS %in% alp.list[[5]]$CHROM_POS) %>%
  filter(CHROM_POS %in% alp.list[[6]]$CHROM_POS) %>%
  filter(CHROM_POS %in% alp.list[[7]]$CHROM_POS) %>%
  filter(CHROM_POS %in% alp.list[[8]]$CHROM_POS)

#making sure all snps are in all pops
#all(snps_alp60md$CHROM_POS %in% alp.list[[6]]$CHROM_POS)

##find highest delta freq abs
hl_alp60md <- separate(snps_alp60md, CHROM_POS, c("SEQ","CHROM", "POS"), sep="_")
hl_alp60md <-hl_alp60md %>% unite('CHROM',SEQ:CHROM, remove = TRUE)
hl_alp60md <- hl_alp60md %>% group_by(CHROM) %>% slice(which.max(DELTA_FREQ_ABS))
hl_alp60md <- hl_alp60md %>% unite('CHROM_POS',CHROM:POS, remove = TRUE)

#final list of cline snps:
clineSNPs_alp60md <- hl_alp60md$CHROM_POS
fwrite(list(clineSNPs_alp60md), file="../clineSNPs_alp60md.txt")

##filter all pops by highest snp (present in all pops) per locus.
alp.list <- lapply(alp.list, function(x) x %>% filter(CHROM_POS %in% clineSNPs_alp60md))

alp_frqs_md60 <- data.frame(alp.list[[1]]$MAJ_FREQ, alp.list[[2]]$MAJ_FREQ, alp.list[[3]]$MAJ_FREQ, alp.list[[4]]$MAJ_FREQ, alp.list[[5]]$MAJ_FREQ, alp.list[[6]]$MAJ_FREQ, alp.list[[7]]$MAJ_FREQ, alp.list[[8]]$MAJ_FREQ)
alp_frqs_md60 <- as.data.frame(t(alp_frqs_md60))
colnames(alp_frqs_md60)<-alp.list[[1]]$CHROM_POS
rownames(alp_frqs_md60)<-NULL

alpsamsis_md60 <- data.frame(alp.list[[1]]$N_CHR, alp.list[[2]]$N_CHR, alp.list[[3]]$N_CHR, alp.list[[4]]$N_CHR, alp.list[[5]]$N_CHR, alp.list[[6]]$N_CHR, alp.list[[7]]$N_CHR, alp.list[[8]]$N_CHR)
alpsamsis_md60 <- as.data.frame(t(alpsamsis_md60))
colnames(alpsamsis_md60)<-paste("N_",alp.list[[1]]$CHROM_POS,sep = "")
rownames(alpsamsis_md60)<-NULL

##alternate allelefrequency and allele sample size columns

alp_frqs_md60 <-cbind.data.frame(alp_frqs_md60,alpsamsis_md60)
alp_frqs_md60 <- as.data.frame(alp_frqs_md60[, c(matrix(1:ncol(alp_frqs_md60), nrow = 2, byrow = T))])

alp_frqs_md60$distance<-alpdist
alp_frqs_md60$population<-popsalp
alp_frqs_md60 <- alp_frqs_md60 %>%
  relocate(distance) %>%
  relocate(population)
write.csv(alp_frqs_md60,"../alp_frqs_md60.csv", row.names = FALSE)

## lapply applies functions to every element of a list. 
pyr.list <- lapply(pyr.list, function(x) x %>% 
                     unite('CHROM_POS',CHROM:POS, remove = TRUE) %>% 
                     separate(MAJ_ALLELE_FREQ, c("MAJ_ALLELE","MAJ_FREQ"), sep=":", remove = TRUE) %>% 
                     separate(MIN_ALLELE_FREQ, c("MIN_ALLELE","MIN_FREQ"), sep=":", remove = TRUE) %>% 
                     filter(CHROM_POS %in% haf_pyr60md_l_onlypos$CHROM_POS))

length(pyr.list[[4]]$CHROM_POS)

snps_pyr60md <- haf_pyr60md_l_onlypos %>% 
  filter(CHROM_POS %in% pyr.list[[1]]$CHROM_POS) %>%
  filter(CHROM_POS %in% pyr.list[[2]]$CHROM_POS) %>%
  filter(CHROM_POS %in% pyr.list[[3]]$CHROM_POS) %>%
  filter(CHROM_POS %in% pyr.list[[4]]$CHROM_POS) %>%
  filter(CHROM_POS %in% pyr.list[[5]]$CHROM_POS) %>%
  filter(CHROM_POS %in% pyr.list[[6]]$CHROM_POS) %>%
  filter(CHROM_POS %in% pyr.list[[7]]$CHROM_POS) %>%
  filter(CHROM_POS %in% pyr.list[[8]]$CHROM_POS) %>%
  filter(CHROM_POS %in% pyr.list[[9]]$CHROM_POS) %>%
  filter(CHROM_POS %in% pyr.list[[10]]$CHROM_POS) %>%
  filter(CHROM_POS %in% pyr.list[[11]]$CHROM_POS)

#making sure all snps are in all pops
#all(snps_pyr60md$CHROM_POS %in% pyr.list[[6]]$CHROM_POS)

##find highest delta freq abs
hl_pyr60md <- separate(snps_pyr60md, CHROM_POS, c("SEQ","CHROM", "POS"), sep="_")
hl_pyr60md <-hl_pyr60md %>% unite('CHROM',SEQ:CHROM, remove = TRUE)
hl_pyr60md <- hl_pyr60md %>% group_by(CHROM) %>% slice(which.max(DELTA_FREQ_ABS))
hl_pyr60md <- hl_pyr60md %>% unite('CHROM_POS',CHROM:POS, remove = TRUE)

#final list of cline snps:
clineSNPs_pyr60md <- hl_pyr60md$CHROM_POS
fwrite(list(clineSNPs_pyr60md), file="../clineSNPs_pyr60md.txt")

##filter all pops by highest snp (present in all pops) per locus.
pyr.list <- lapply(pyr.list, function(x) x %>% filter(CHROM_POS %in% clineSNPs_pyr60md))

pyr_frqs_md60 <- data.frame(pyr.list[[1]]$MAJ_FREQ, pyr.list[[2]]$MAJ_FREQ, pyr.list[[3]]$MAJ_FREQ, pyr.list[[4]]$MAJ_FREQ, pyr.list[[5]]$MAJ_FREQ, pyr.list[[6]]$MAJ_FREQ, pyr.list[[7]]$MAJ_FREQ, pyr.list[[8]]$MAJ_FREQ, pyr.list[[9]]$MAJ_FREQ, pyr.list[[10]]$MAJ_FREQ, pyr.list[[11]]$MAJ_FREQ)
pyr_frqs_md60 <- as.data.frame(t(pyr_frqs_md60))
colnames(pyr_frqs_md60)<-pyr.list[[1]]$CHROM_POS
rownames(pyr_frqs_md60)<-NULL

pyrsamsis_md60 <- data.frame(pyr.list[[1]]$N_CHR, pyr.list[[2]]$N_CHR, pyr.list[[3]]$N_CHR, pyr.list[[4]]$N_CHR, pyr.list[[5]]$N_CHR, pyr.list[[6]]$N_CHR, pyr.list[[7]]$N_CHR, pyr.list[[8]]$N_CHR, pyr.list[[9]]$N_CHR, pyr.list[[10]]$N_CHR, pyr.list[[11]]$N_CHR)
pyrsamsis_md60 <- as.data.frame(t(pyrsamsis_md60))
colnames(pyrsamsis_md60)<-paste("N_",pyr.list[[1]]$CHROM_POS,sep = "")
rownames(pyrsamsis_md60)<-NULL

##alternate allelefrequency and allele sample size columns

pyr_frqs_md60 <-cbind.data.frame(pyr_frqs_md60,pyrsamsis_md60)
pyr_frqs_md60 <- as.data.frame(pyr_frqs_md60[, c(matrix(1:ncol(pyr_frqs_md60), nrow = 2, byrow = T))])

pyr_frqs_md60$distance<-pyrdist
pyr_frqs_md60$population<-popspyr
pyr_frqs_md60 <- pyr_frqs_md60 %>%
  relocate(distance) %>%
  relocate(population)
write.csv(pyr_frqs_md60,"../pyr_frqs_md60.csv", row.names = FALSE)

############################################
############################################
########   Convert into HZAR input  ########
############################################
############################################

pyr60mdmodelSNP1FitR <-
  hzar.first.fitRequest.old.ML(model=pyr_md60.models[[1]] ,
                               pyr_md60.SNPlist[[1]],
                               verbose=TRUE);
pyr60mdmodelSNP1FitR$mcmcParam$chainLength <- 2e3;
pyr60mdmodelSNP1FitR$mcmcParam$burnin <- 5e2;
pyr60mdmodelSNP1Fit <- hzar.doFit(pyr60mdmodelSNP1FitR)
plot(hzar.mcmc.bindLL(pyr60mdmodelSNP1Fit))
pyr60mdmodelSNP1Data <-
  hzar.dataGroup.add(pyr60mdmodelSNP1Fit);
## Not run: 
pyr60mdmodelSNP1Data <-
  hzar.dataGroup.add(
    pyr60mdmodelSNP1Data,
    hzar.chain.doSeq(hzar.next.fitRequest(pyr60mdmodelSNP1Fit)));
hzar.plot.cline(pyr60mdmodelSNP1Data);
hzar.plot.fzCline(pyr60mdmodelSNP1Data);

## End(Not run)
print(hzar.getLLCutParam(pyr60mdmodelSNP1Data,c("center","width")));
pyr60mdmodelSNP1Null <- hzar.dataGroup.null(pyr_frqs_md60_SNP1);
pyr_frqs_md60_SNP1dGs <- list(clineModel = pyr60mdmodelSNP1Data,
                   nullModel  = pyr60mdmodelSNP1Null);
pyr_frqs_md60_SNP1oDG <- hzar.make.obsDataGroup(pyr_frqs_md60_SNP1dGs);
pyr_frqs_md60_SNP1oDG <- hzar.copyModelLabels(pyr_frqs_md60_SNP1dGs,pyr_frqs_md60_SNP1oDG);
hzar.plot.cline(pyr_frqs_md60_SNP1oDG);
print(hzar.AICc.hzar.obsDataGroup(pyr_frqs_md60_SNP1oDG));

















### read in prefiltered vcf file from

head(vcf_alp_md40_l)

### get stats: 
vcf_alp_md40_l
# ***** Object of Class vcfR *****
#144 samples
#485 CHROMs
#9,828 variants
#Object size: 19.5 Mb
#0 percent missing data
# *****        *****         *****


#must set parallel.core = 1 or error: 
#One of the nodes produced an error: Can not open file .gds'. 
#The process cannot access the file because it is being used by another process.
alp40VCF <- tidy_vcf(data="alp_md40_l.recode.vcf", parallel.core = 1)
hzar_filtered <- write_hzar(
  tidyVCF,
  distances = "distances2.txt",
  filename = NULL,
  parallel.core = parallel::detectCores() - 1
)


genomic_converter("alp_md40_l.recode.vcf", parallel.core = 1, output = "hzar_alp_md40_l.hzar", verbose = TRUE)
genomic_converter("alp_md60_l.recode.vcf",output = "hzar_alp_md60_l.hzar", verbose = TRUE)
genomic_converter("pyr_md40_l.recode.vcf",output = "hzar_pyr_md40_l.hzar", verbose = TRUE)
genomic_converter("pyr_md60_l.recode.vcf",output = "hzar_pyr_md60_l.hzar", verbose = TRUE)

write_hzar(conv_data,"hzar_alp_m40_l.hzar")
write_hzar(conv_data,"hzar_alp_m60_l.hzar")
write_hzar(conv_data,"hzar_pyr_m40_l.hzar")
write_hzar(conv_data,"hzar_pyr_m60_l.hzar")