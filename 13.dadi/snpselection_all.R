### make useful funcion
sampleWithoutSurprises <- function(x) {
  if (length(x) <= 1) {
    return(x)
  } else {
    return(sample(x,1))
  }
}

### set seed for reproducibility
set.seed(6969)

### set working directory
getwd()
setwd(".")

colns <-c("CHROM","POS")

### read broken bcf file, skip header (47 rows), ignore variant calls and other rows
snpsfile <- read.csv("01.output/dadi_all.bcf", sep="\t", header=FALSE, row.names = NULL, skip = 47, colClasses=c(NA, NA, rep("NULL",118)))

### get linked snps from this
write(snpsfile$V2,"linkedsnps_all.txt")

### read list of all loci
locilist<-read.csv("../00.input/bait.info", sep="\t", header=FALSE, row.names = NULL)

snps<-list()

### make list of random snps selected from locus (only if locus has snps)
for (i in 1:length(locilist$V2)){
  if (!all(snpsfile$V2 %in% locilist[[i,2]]:locilist[[i,3]]==FALSE)){
    k<-sampleWithoutSurprises(which(snpsfile$V2 %in% locilist[[i,2]]:locilist[[i,3]]));
    snps[[i]]<-snpsfile[[k,2]]
  }
}
snps<-unlist(snps, use.names = FALSE)
write(snps,"randomsnps_all.txt")
