### load packages
library(tidyverse)
library(ggplot2)
library(data.table)
library(hzar)
library(dplyr)

### set to working directory
getwd()
setwd("C:/Users/Jag/Documents/LMU/work/01.capture/07.hzar/pop_mafs")

### load allele freq files in "/pop_mafs/"
dists<- c(0,5576,8853,12168,13494,15183,17845,18700,20818,23140,26761,29176)
pops<- c("GAB","HER", "SOQ", "TOU", "ARA", "POR", "MUL", "TRO", "FOR", "PAZ","LAN", "ESC")
colns <-c("CHROM","POS","N_ALLELES","N_CHR","MAJ_ALLELE_FREQ","MIN_ALLELE_FREQ") 


af.list<-list(read.csv("vcf_af_GAB.frq", sep="\t", header=TRUE, row.names =  NULL, col.names = colns),
              read.csv("vcf_af_HER.frq", sep="\t", header=TRUE, row.names =  NULL, col.names = colns),
              read.csv("vcf_af_SOQ.frq", sep="\t", header=TRUE, row.names =  NULL, col.names = colns),
              read.csv("vcf_af_TOU.frq", sep="\t", header=TRUE, row.names =  NULL, col.names = colns),
              read.csv("vcf_af_ARA.frq", sep="\t", header=TRUE, row.names =  NULL, col.names = colns),
              read.csv("vcf_af_POR.frq", sep="\t", header=TRUE, row.names =  NULL, col.names = colns),
              read.csv("vcf_af_MUL.frq", sep="\t", header=TRUE, row.names =  NULL, col.names = colns),
              read.csv("vcf_af_TRO.frq", sep="\t", header=TRUE, row.names =  NULL, col.names = colns),
              read.csv("vcf_af_FOR.frq", sep="\t", header=TRUE, row.names =  NULL, col.names = colns),
              read.csv("vcf_af_PAZ.frq", sep="\t", header=TRUE, row.names =  NULL, col.names = colns),
              read.csv("vcf_af_LAN.frq", sep="\t", header=TRUE, row.names =  NULL, col.names = colns),
              read.csv("vcf_af_ESC.frq", sep="\t", header=TRUE, row.names =  NULL, col.names = colns))

af.list<-lapply(af.list, function(x) x %>% 
                  separate(MAJ_ALLELE_FREQ, c("MAJ_ALLELE","MAJ_FREQ"), sep=":", remove = TRUE) %>% 
                  separate(MIN_ALLELE_FREQ, c("MIN_ALLELE","MIN_FREQ"), sep=":", remove = TRUE) %>% 
                  dplyr::select(POS,N_CHR,MAJ_ALLELE,MAJ_FREQ,MIN_ALLELE,MIN_FREQ))

## major allele is the same across all pops :)

af_frqs <- data.frame(af.list[[1]]$MAJ_FREQ, af.list[[2]]$MAJ_FREQ, af.list[[3]]$MAJ_FREQ, af.list[[4]]$MAJ_FREQ, af.list[[5]]$MAJ_FREQ, af.list[[6]]$MAJ_FREQ, af.list[[7]]$MAJ_FREQ, af.list[[8]]$MAJ_FREQ, af.list[[9]]$MAJ_FREQ, af.list[[10]]$MAJ_FREQ, af.list[[11]]$MAJ_FREQ, af.list[[12]]$MAJ_FREQ)
af_frqs <- as.data.frame(t(af_frqs))
colnames(af_frqs)<-af.list[[1]]$POS
rownames(af_frqs)<-NULL

af_nchrs <- data.frame(af.list[[1]]$N_CHR, af.list[[2]]$N_CHR, af.list[[3]]$N_CHR, af.list[[4]]$N_CHR, af.list[[5]]$N_CHR, af.list[[6]]$N_CHR, af.list[[7]]$N_CHR, af.list[[8]]$N_CHR, af.list[[9]]$N_CHR, af.list[[10]]$N_CHR, af.list[[11]]$N_CHR, af.list[[12]]$N_CHR)
af_nchrs <- as.data.frame(t(af_nchrs))
colnames(af_nchrs)<-paste("N_",af.list[[1]]$POS,sep = "")
rownames(af_nchrs)<-NULL

##alternate allelefrequency and allele sample size columns

af_frqs <-cbind.data.frame(af_frqs,af_nchrs)
af_frqs <- as.data.frame(af_frqs[, c(matrix(1:ncol(af_frqs), nrow = 2, byrow = T))])

af_frqs$distance<-dists
af_frqs$population<-pops
af_frqs <- af_frqs %>%
  relocate(distance) %>%
  relocate(population)
write.csv(af_frqs,"../df_frqs.csv", row.names = FALSE)

rm(af.list)
############################################
############################################
##########  HZAR CLINE MODELLING  ##########
############################################
############################################

## A typical chain length.  This value is the default setting in the package.
chainLength=2e5;
burnin=5e4;

## Make each model run off a separate seed
mainSeed=
  list(A=c(596,528,124,978,544,99),
       B=c(528,124,978,544,99,596),
       C=c(124,978,544,99,596,528))


### load hzar formatted snp frequency dataset
af_frqs <- read.csv("../df_frqs.csv")

af_frqs.SNPs<- list()

## chose columns that contain snp frequencies
snplist<- seq(3,ncol(af_frqs),2)

# polarise: ascending

for (i in snplist){
  if (af_frqs[1,i] > af_frqs[11,i]){
    af_frqs[[i]]<- sapply(af_frqs[[i]], function(x) 1-x)
  }
}


## make a list of the snp data entered into hzar, each element containing the snp, the distance vector, the snp information. max 15 min (8000 snps)
for (i in snplist){
  #print(i)
  k<-((i-1)/2)
  j<-i+1
  SNP<- paste0("SNP",k)
  tmp<- af_frqs %>% 
    dplyr::select(c(2,all_of(i),all_of(j)))
  tmp<- hzar.doMolecularData1DPops(as.numeric(unlist(tmp[1])),as.numeric(unlist(tmp[2])),as.numeric(unlist(tmp[3])))
  af_frqs.SNPs[[SNP]]$obs<-tmp
}
##remove loop residue
remove(tmp,i,j,k,SNP,af_frqs,af_nchrs,snplist)

af_frqs.SNPs<- lapply(af_frqs.SNPs, function(x){
  x$models<-list();
  x$fitR<-list();
  x$runs<-list();
  x$analysis<-list();
  return(x)
})

## define models: easiest sigmoidal function and define limits of parameter space
for (i in 1:length(af_frqs.SNPs)){
  af_frqs.SNPs[[i]]$models<-hzar.makeCline1DFreq(af_frqs.SNPs[[i]]$obs);
  af_frqs.SNPs[[i]]$models<-hzar.model.addBoxReq(af_frqs.SNPs[[i]]$models,low = -5,high = 30005)
}

saveRDS(af_frqs.SNPs, file="af_frqs_SNPs")
af_frqs.SNPs<-readRDS("af_frqs_SNPs")

## fit request the model to each snp (4h)
for (i in 1:length(af_frqs.SNPs)){
  af_frqs.SNPs[[i]]$fitR$init<-hzar.first.fitRequest.old.ML(model = af_frqs.SNPs[[i]]$models,obsData = af_frqs.SNPs[[i]]$obs, verbose = FALSE)
}

af_frqs.SNPs<-lapply(
  af_frqs.SNPs,function(x){
    x$fitR$init$mcmcParam$chainLength <-chainLength; #2e5
    x$fitR$init$mcmcParam$burnin <-chainLength %/% 10;
    x$fitR$init$mcmcParam$seed[[1]] <-mainSeed$A;
    return(x)
  }
)

#pyr_md40md.models.FitR <- mapply(function(x,y){
#  lapply(1:354, function(el) hzar.first.fitRequest.old.ML(model=x ,obsData = y, verbose=TRUE))
#},x=pyr_md40.models ,y=pyr_md40.SNPlist)

saveRDS(af_frqs.SNPs, file="af_frqs_SNPs")
af_frqs.SNPs<-readRDS("af_frqs_SNPs")

## run initital chains to compile new fitrequest (14:25 - 14:58 --> 24:05 - 13:30)

af_frqs.SNPs <- lapply(af_frqs.SNPs, function(x){ 
  x$runs$init <- list();
  x$runs$init <- hzar.doFit(x$fitR$init)
  return(x)
})


## Compile a new set of fit requests using the initial chains (13:45 - ? max 40 min)
af_frqs.SNPs <- lapply(af_frqs.SNPs, function(x){ x$fitR$chains<-hzar.next.fitRequest(x$runs$init);return(x)})

saveRDS(af_frqs.SNPs, file="prep_frqs_SNPs")
af_frqs.SNPs <- readRDS("prep_frqs_SNPs")


## Replicate each fit request 3 times, keeping the original
## seeds while switching to a new seed channel.
af_frqs.SNPs <- lapply(af_frqs.SNPs, function (x){ x$fitR$chains<-hzar.multiFitRequest(x$fitR$chains, each=3, baseSeed=NULL); return(x)})

## runif(3,-5,30005) center
af_frqs.SNPs <- lapply(af_frqs.SNPs, function (x){ 
  x$fitR$chains[[1]]$modelParam$init["center"]= 10244.295;
  x$fitR$chains[[2]]$modelParam$init["center"]= 14244.592;
  x$fitR$chains[[3]]$modelParam$init["center"]= 3543.751;
  
  ## runif(3,0,30010) width 
  x$fitR$chains[[1]]$modelParam$init["width"]= 9573.702; 
  x$fitR$chains[[2]]$modelParam$init["width"]= 3514.500;  
  x$fitR$chains[[3]]$modelParam$init["width"]= 24517.841;
  return(x)
})

## Go ahead and run a chain of 3 runs for every fit request
#af_frqs.SNPs <-  lapply(af_frqs.SNPs,function(x){
#  x$runs$chains <- hzar.doChain.multi(x$fitR$chains,doPar=TRUE,inOrder=FALSE,count=3);
#  return(x)
#})

firstpart_frqs.SNPs<-af_frqs.SNPs[1:3006]
secondpart_frqs.SNPs<-af_frqs.SNPs[3007:6000]
thirdpart_frqs.SNPs<-af_frqs.SNPs[6001:length(af_frqs.SNPs)]

firstpart_frqs.SNPs <-  lapply(firstpart_frqs.SNPs,function(x){
  x$runs$chains <- hzar.doChain.multi(x$fitR$chains,doPar=TRUE,inOrder=FALSE,count=3);
  return(x)
})

saveRDS(firstpart_frqs.SNPs,file="firstpart_frqs_SNPs")
rm(firstpart_frqs.SNPs)

secondpart_frqs.SNPs <-  lapply(secondpart_frqs.SNPs,function(x){
  x$runs$chains <- hzar.doChain.multi(x$fitR$chains,doPar=TRUE,inOrder=FALSE,count=3);
  return(x)
})

saveRDS(secondpart_frqs.SNPs,file="secondpart_frqs_SNPs")
rm(secondpart_frqs.SNPs)

thirdpart_frqs.SNPs <-  lapply(thirdpart_frqs.SNPs,function(x){
  x$runs$chains <- hzar.doChain.multi(x$fitR$chains,doPar=TRUE,inOrder=FALSE,count=3);
  return(x)
})

saveRDS(thirdpart_frqs.SNPs,file="thirdpart_frqs_SNPs")
rm(thirdpart_frqs.SNPs)

firstpart_frqs.SNPs<-readRDS("firstpart_frqs_SNPs")
secondpart_frqs.SNPs<-readRDS("secondpart_frqs_SNPs")
thirdpart_frqs.SNPs<-readRDS("thirdpart_frqs_SNPs")

all_frqs.SNPs<-c(firstpart_frqs.SNPs,secondpart_frqs.SNPs,thirdpart_frqs.SNPs)
saveRDS(all_frqs.SNPs,file="all_frqs_SNPs")


rm(firstpart_frqs.SNPs,secondpart_frqs.SNPs,thirdpart_frqs.SNPs)
#did it converge?
all_frqs.summary<-list()


for (i in 1:length(all_frqs.SNPs)){
  all_frqs.summary[[i]]<-summary(do.call(mcmc.list,
                                         lapply(all_frqs.SNPs[[i]]$runs$chains[1:3],
                                                function(x) hzar.mcmc.bindLL(x[[3]]) )) )
}

saveRDS(all_frqs.summary, file="all_frqs_summary")
all_frqs.summary <- readRDS("all_frqs_summary")

rm(all_frqs.SNPs)
# extract centre values
all_frqs.centres<-list()
for (i in 1:length(all_frqs.summary)){
  all_frqs.centres[[i]]<-all_frqs.summary[[i]]$statistics[[1]]
}

all_frqs.width<-list()
for (i in 1:length(all_frqs.summary)){
  all_frqs.width[[i]]<-all_frqs.summary[[i]]$statistics[[2]]
}

cent_widths<-as.data.frame(unlist(all_frqs.centres))
cent_widths$width<-(unlist(all_frqs.width))
colnames(cent_widths)<-"centre"


ggplot(cent_widths, aes(x=centre)) +
  geom_density(fill="#69b3a2", color="#e9ecef", alpha=0.6) +
  theme_classic()

ggplot(cent_widths, aes(x=width)) +
  geom_density(fill="#69b3a2", color="#e9ecef", alpha=0.6) +
  theme_classic()


ggplot(cent_widths, aes(centre))+ 
  geom_histogram(binwidth=500, fill="#69b3a2") +
  theme_bw()

ggplot(cent_widths, aes(width))+ 
  geom_histogram(binwidth=500, fill="#69b3a2") +
  theme_bw()

ggplot(cent_widths, aes(width,centre)) +
  geom_point()

qplot(as.factor(cent_widths), geom="histogram")

## Start aggregation of data for analysis

## Create a model data group for the null model (expected allele
## frequency independent of distance along cline) to include in
## analysis.
for (i in 1:length(pyr_md40.SNPs)){
  pyr_md40.SNPs[[i]]$analysis$initDGs <- list(
    nullModel =  hzar.dataGroup.null(pyr_md40.SNPs[[1]]$obs))
}

## Create a model data group (hzar.dataGroup object) for each
## model from the initial runs.
for (i in 1:length(pyr_md40.SNPs)){
  pyr_md40.SNPs[[i]]$analysis$initDGs$model1 <- hzar.dataGroup.add(pyr_md40.SNPs[[i]]$runs$init)
}

## Create a hzar.obsDataGroup object from the four hzar.dataGroup
## just created, copying the naming scheme (nullModel, modelI).
## skip, because not possible for multiple loci. 
## "Error in hzar.make.obsDataGroup(otherDataGroups, obsDataGroup) : 
## All dataGroups must be from the same observation data"
## no workaround found. 

for (i in 1:length(pyr_md40.SNPs)){ 
  pyr_md40.SNPs[[i]]$analysis$oDG <- hzar.make.obsDataGroup(pyr_md40.SNPs[[i]]$analysis$initDGs);
  pyr_md40.SNPs[[i]]$analysis$oDG <- hzar.copyModelLabels(pyr_md40.SNPs[[i]]$analysis$initDGs,pyr_md40.SNPs[[i]]$analysis$oDG)
}

pyr_md40.SNPs[[1]]$analysis$oDG <- hzar.make.obsDataGroup(lapply(pyr_md40.SNPs[[1]]$runs$chains, 
                                                                 hzar.dataGroup.add),
                                                          pyr_md40.SNPs[[1]]$analysis$oDG)

pyr_md40.SNPs[[2]]$analysis$oDG <- hzar.make.obsDataGroup(pyr_md40.SNPs[[2]]$analysis$initDGs)

mkn$AdaA$analysis$oDG <-
  hzar.make.obsDataGroup(mkn$AdaA$analysis$initDGs)
mkn$AdaA$analysis$oDG <-
  hzar.copyModelLabels(mkn$AdaA$analysis$initDGs,
                       mkn$AdaA$analysis$oDG)

## Convert all 27 runs to hzar.dataGroup objects, adding them to
## the hzar.obsDataGroup object.
pyr_md40.SNPs[[1]]$analysis$oDG <-
  hzar.make.obsDataGroup(lapply(pyr_md40.SNPs[[1]]$runs$chains,
                                hzar.dataGroup.add),
                         pyr_md40.SNPs[[1]]$analysis$oDG);