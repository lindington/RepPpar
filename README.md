# target/capture including the new samples I sent for sequencing in dec 2021 -> data from 2022

I contributed to the following components of this project: 
- [Assembly](/01.paleomix/)
- [PCA](/04.pca)
- [NGSadmix](/05.ngsadmix/)
- [Fst](/08.fst/)
- [Heterozygosity](/09.hz/)
- [Pi, Theta, and Tajima's D](/10.pi_theta_taj/)
- [Geographic Clines](/11.hzar/)
- [Dxy](/12.dxy/)
- [Demographic Analysis](/13.dadi/)
- [Genomic clines](/14.gghybrid)

> **NOTE**: MUL414 has to be reassigned to SOQ414 throughout the analysis! (shouldn't rename bamfiles)

## Assembly
#### Paleomix

I used the [BAM Pipeline](https://paleomix.readthedocs.io/en/stable/bam_pipeline/index.html) from Paleomix (Schubert et al., 2014) to create ``.bam`` files (compressed binary version of a SAM file, i.e. text-based format originally for storing biological sequences aligned to a reference sequence), on the basis of enriques scripts. See the relevant files in [01.paleomix](/01.paleomix/README.md).

#### Depth

I calculated depth and coverage. See the relevant files in [02.depth](02.depth/).

## PCA and NGSadmix

I started doing population structure analyses using the full dataset including a geographical "outgroup" (bamlist_out.txt). To do any analysis in angsd, I needed a bamlist (`.txt` files containing paths to relevant bamfiles), which i made using `readlink -f ../01.paleomix/*bam > ../00.input/bamlist122.txt`. 

I started by making a `.geno` and a `.beagle` file with the same filters. The ``.beagle`` file will be used for [NGSadmix](#ngsadmix). I specified to retain the baited regions only used the `-sites` flag and restricted analyses to chr1 using the `-r` flag.

I ran the following code on the module installation of angsd (v. 0.933-gcc8)

```bash
#! /bin/bash

#SBATCH -J beagle
#SBATCH --output=beagle.out
#SBATCH --clusters=biohpc_gen
#SBATCH --partition=biohpc_gen_production
#SBATCH --cpus-per-task=40
#SBATCH --mem=600gb
#SBATCH --mail-type=ALL
#SBATCH --mail-user=hagberg@lmu.de
#SBATCH -t 4-00:00:00

STARTTIME=$(date +"%s")
echo " Job started at $(date)."

module load angsd/0.933-gcc8

angsd -b /dss/dsslegfs01/pr53da/pr53da-dss-0029/capture/00.input/bamlist_out.txt -ref /dss/dsslegfs01/pr53da/pr53da-dss-0029/TranscriptomeAnalyses/TranscriptomeReference/grasshopperWolbRef.fasta -doMajorMinor 1 -GL 1 -doGlf 2 -SNP_pval 1e-2 -doMaf 1 -nThreads 4 -r chr1: -sites /dss/dsslegfs01/pr53da/pr53da-dss-0029/capture/00.input/bait.sites -baq 1 -remove_bads 1 -uniqueOnly 1 -C 50 -minMapQ 15 -only_proper_pairs 0 -minQ 20 -doCounts 1 -doPost 2 -doGeno 32 -minInd 176 -setMaxDepth 14472 -out 01.output/220_out_maxdepth

ENDTIME=$(date +%s)
TIMESPEND=$(($ENDTIME - $STARTTIME))
((sec=TIMESPEND%60, TIMESPEND/=60, min=TIMESPEND%60, hrs=TIMESPEND/60))
timestamp=$(printf "%d:%02d:%02d" $hrs $min $sec)
echo "Job ended at $(date). Took $timestamp hours:minutes:seconds to complete."
```

#### PCA 

To make a PCA, I unzipped the `.geno.gz` file using `gzip -d .geno` and counted the number of variable sites using `zcat ../beagle/01.output/beagle_bait.mafs.gz | tail -n+2 | wc -l`. There are 371391 variable sites.

I then used ngsCovar to make a covariance table. To use ngsCovar, you have to intall `ngsTools` LOCALLY. I will ask IT to install it as a module. There is no conda package for it, and ngsCovar is not included in ANGSD. In the past, I've had issues with the dependencies being incompatible (htslib vs ngstools or something). This was not an issue this time.

To make the covariance table, i just ran the following command on bash.
```bash
/dss/dsshome1/lxc0A/ru67vil/programs/ngsTools/ngsPopGen/ngsCovar -probfile ../beagle/01.output/beagle_bait.geno -outfile beagle_bait.covar -nind 98 -nsites 371391 -call 0 -norm 0` 
```
Then, I made a plink cluster file using 

```bash
Rscript -e 'write.table(cbind(seq(1,50),rep(1,50),c(rep("ARA",10),rep("ESC",5),rep("FOR",8),rep("GAB",5),rep("HER",10),rep("LAN",10),rep("MUL",7),rep("SOQ",1),rep("MUL",2),rep("PAZ",10),rep("POR",10),rep("SOQ",10),rep("TOU",10))),row.names=F,col.names=c("FID","IID","CLUSTER"),sep ="\t",file="plink_bait98.clst",quote=F)'
```
These assignments have to follow the order of the ``.bam`` files in the bamlist used to make the ``.beagle`` file.
Note that nine MUL assignments are separated by a single assignment to SOQ. This is because the bamfile of a mislabeled individual from Soques was still grouped with Mulas individuals. 

To plot the PCA components, I ran the following script locally on my laptop, because the R packages needed aren't installed on the server. 

```R
library(ggplot2)

# Annotation file is in plink cluster format
# Read input file
setwd("C:/Users/Jag/Documents/LMU/work/01.capture/00.baitanalyses/")
covar <- read.table('beagle_bait.covar', stringsAsFactors = FALSE);
# Read annot file
annot <- read.table('plink_bait98.clst',sep="\t",header=TRUE);
#note that plink cluster files are usually tab - separated instead
#View(annot)
# Parse components to analyze
comp <- as.numeric(strsplit('1-2',"-",fixed=TRUE)[[1]])
View(comp[1])
# Eigenvalues
eig <- eigen(covar,symm=TRUE);
eig$val <- eig$val/sum(eig$val);
cat(signif(eig$val,digits=3)*100,"\n");

# Plot
PC <- as.data.frame(eig$vectors)
colnames(PC) <- gsub("V","PC",colnames(PC))
PC$Pop <- factor(annot$CLUSTER)
PC$Col <- factor(annot$COLOUR)
#View(PC)

title <- paste("PC",comp[1],"(",signif(eig$val[comp[1]],digits=3)*100,"%)","/PC",comp[2],"(",signif(eig$val[comp[2]],digits=3)*100,"%)",sep="",collapse="")
#View(title)
x_axis = paste("PC",comp[1],sep ="")
y_axis = paste("PC",comp[2],sep ="")

#colours set to match throughout the figures of the paper
#colperpop <- c("GAB"="#FFCC00", "HER"="#CCFF66", "SOQ"="#CCFFCC", "TOU"="#66FF99", "ARA"="#99FFFF", "POR"="#6699FF", "MUL"="#99CCFF", "FOR"="#FFCCFF", "PAZ"="#FF99FF", "LAN"="#FF3399", "ESC"="#FF3333")

colperpop <- c("GAB"="#001219", "HER"="#305762", "SOQ"="#447B8A", "TOU"="#79968E", "ARA"="#BEC6A5", "POR"="#FFF5BE", "MUL"="#F2D18A", "FOR"="#E4AD5B", "PAZ"="#CF7B4E", "LAN"="#C23436", "ESC"="#8D1216")
#colperpop <- c("ARA"="#000000","ESC"="#5C4B51","FOR"="#D10000","GAB"="#8CBEB2","HER"="#F2EBBF","LAN"="#F3B562","MUL"="#63a02c","PAZ"="#589CDA","POR"="#A4A4A4","SOQ"="#26549C","TOU"="#489CDA")


ggplot() + 
  geom_point(aes_string(x=x_axis,y=y_axis, fill=PC$Pop),
             data=PC,
             colour="white",
             pch = 21,  
             size = 3, 
             stroke = 0.4, 
             show.legend = TRUE) +
  scale_y_continuous(minor_breaks=NULL, limits=c(-0.2,0.2)) +
  scale_x_continuous(minor_breaks=NULL, limits=c(-0.2,0.2)) + 
  coord_fixed()  +
  ggtitle(title) + 
  scale_fill_manual(values=colperpop,name="Populations",breaks=c("GAB","HER","SOQ","TOU","ARA","POR","MUL","FOR","PAZ","LAN","ESC"),labels=c("Gabas","Hermine","Soques","Tourmont","Araille","Portalet","Mulas","Formigal","Pazino","Lanuza","Escarilla"))+
  
  ggsave('bait_PCA12.pdf')

unlink("Rplots.pdf",force=TRUE)

###################################
#### Start components 3 and 4: ####
###################################
comp <- as.numeric(strsplit('3-4',"-",fixed=TRUE)[[1]])
View(comp[1])

eig <- eigen(covar,symm=TRUE);
eig$val <- eig$val/sum(eig$val);
cat(signif(eig$val,digits=3)*100,"\n");

# Plot
PC <- as.data.frame(eig$vectors)
colnames(PC) <- gsub("V","PC",colnames(PC))
PC$Pop <- factor(annot$CLUSTER)
PC$Col <- factor(annot$COLOUR)
#View(PC)


title <- paste("PC",comp[1],"(",signif(eig$val[comp[1]],digits=3)*100,"%)","/PC",comp[2],"(",signif(eig$val[comp[2]],digits=3)*100,"%)",sep="",collapse="")
View(title)
x_axis = paste("PC",comp[1],sep ="")
y_axis = paste("PC",comp[2],sep ="")


ggplot() + 
  geom_point(aes_string(x=x_axis,y=y_axis, fill=PC$Pop),
             data=PC,
             colour="white",
             pch = 21,  
             size = 3, 
             stroke = 0.4, 
             show.legend = TRUE) +
  #scale_y_continuous(minor_breaks=NULL, limits=c(-0.2,0.2)) +
  #scale_x_continuous(minor_breaks=NULL, limits=c(-0.2,0.2)) + 
  #coord_fixed()  +
  ggtitle(title) + 
  scale_fill_manual(values=colperpop,name="Populations",breaks=c("GAB","HER","SOQ","TOU","ARA","POR","MUL","FOR","PAZ","LAN","ESC"),labels=c("Gabas","Hermine","Soques","Tourmont","Araille","Portalet","Mulas","Formigal","Pazino","Lanuza","Escarilla"))+
  
  ggsave('bait_PCA34.pdf')

unlink("Rplots.pdf",force=TRUE)

```

### EMU-PCA
Because of the large amount of missing data in some individuals, I decided to additionally make an EMU-PCA, which imputes missing data better for extremely large amounts, designed especially for target/capture data.
However, EMU-PCA uses plink style ``.bed`` files and cannot be performed on genotype likelihood files (e.g. ``.beagle`` extention gl files). Since there is no pipeline for using ANGSD files for emu-pca input, I tried different workarounds.

1. I used the beta option ``-doBcf 1`` in the same script that i produced the ``.beagle`` file with, to see if i can produce an adequate VCF file for downstream analyses. This then has to be converted to plink input format.

2. Beagle has some utilities i can use in combination to produce a vcf file with this pipeline: 
   - ``gprobs2beagle.jar`` which calls genotypes into a ``.beagle`` file
   - ``beagle2vcf.jar`` which converts the genotype file into a ``.vcf`` file
   - plink has an option to convert ``.vcf`` files into plink style ``.bed`` files.

3. I can check the options for using ANGSD directly to call genotypes using ``-doGeno 2`` (from ANGSD: ``"2: write the called genotype encoded as -1,0,1,2, -1=not called"``)
 
> skipped, too much effort 
> and overkill

### Depth for excluding MUL ind 
#### Using ANGSD
I used the ANGSD option ``-doDepth`` to calculate coverage of individual samples to assess whether any individual had extremely low coverage. All the quality filters were kept the same as when making the ``.beagle`` and ``.geno`` files. I used the following script:

```bash
#!/bin/bash
#SBATCH -J depth_
#SBATCH --output=depthind.out
#SBATCH --clusters=biohpc_gen
#SBATCH --partition=biohpc_gen_production
#SBATCH --cpus-per-task=40
#SBATCH --mem=600gb
#SBATCH --mail-type=ALL
#SBATCH --mail-user=hagberg@lmu.de
#SBATCH -t 4-00:00:00

STARTTIME=$(date +"%s")

module load angsd/0.933-gcc8

angsd -b /dss/dsslegfs01/pr53da/pr53da-dss-0029/capture/00.input/bamlist98.txt -ref /dss/dsslegfs01/pr53da/pr53da-dss-0029/TranscriptomeAnalyses/TranscriptomeReference/grasshopperWolbRef.fasta -doDepth 1 -out 01.output/depth_ind -doCounts 1 -r chr1: -minMapQ 15 -minQ 20 -nThreads 4 -sites /dss/dsslegfs01/pr53da/pr53da-dss-0029/capture/00.input/bait.sites -baq 1 -remove_bads 1 -uniqueOnly 1 -C 50 -only_proper_pairs 0  -minInd 79

ENDTIME=$(date +%s)
TIMESPEND=$(($ENDTIME - $STARTTIME))
((sec=TIMESPEND%60, TIMESPEND/=60, min=TIMESPEND%60, hrs=TIMESPEND/60))

timestamp=$(printf "%d:%02d:%02d" $hrs $min $sec)

echo "Job ended at $(date). Took $timestamp hours:minutes:seconds to complete."
```

The output is given as individual and whole population depth. 
Coverage high in some genes, representing TEs?

>Output with mandatory maximum depth means no accurate representation of highest depth, or distribution of regions with high depth -> problematic because no full picture of highest depth.

- tried different max depths to make sure the last column is always 0 (higher depth gets counted as maxdepth gene)
- tried: 100 (fail, most regions have higher depth)
- tried: 1000000000 (fail, i.e. files unopenable with vim)
- tried: 1000000 (works, files openable with vim, last column always 0).

-> use maxvim = 1000000, but not elegant.

I adapted Vitalis scripts for quick plotting.

>-dodepth calculates perbase and you cant sum over regions. 
> continue with bedtools.

#### Using bedtools

I used bedtools to calculate the global depth per baited region, global depth per base pair, depth per individual, and depth per individual per base pair: with the following script:

```bash
#! /bin/bash

#SBATCH --clusters=biohpc_gen
#SBATCH --partition=biohpc_gen_normal
#SBATCH --mail-user=hagberg@lmu.de
#SBATCH --mail-type=FAIL
#SBATCH --cpus-per-task=10
#SBATCH --mem-per-cpu=2000mb
#SBATCH --time=2-00:00:00
#SBATCH -J GlobBP
#SBATCH -o globbp.out

STARTTIME=$(date +%s)

echo " ###### Job name: $SLURM_JOB_NAME, job ID: $SLURM_JOB_ID ######"
echo " Job started at $(date)."

module load bedtools2/2.27.1

for ind in ARA270 ARA271 ARA273 ARA275 ARA276 ARA277 ARA279 ARA280 ARA281 ARA285 ESC011 ESC012 ESC013 ESC015 ESC344 FOR876 FOR877 FOR881 FOR884 FOR885 FOR887 FOR888 FOR890 GAB512 GAB513 GAB514 GAB515 GAB517 HER450 HER451 HER452 HER453 HER454 HER456 HER457 HER459 HER464 HER465 LAN927 LAN928 LAN929 LAN930 LAN931 LAN934 LAN935 LAN936 LAN938 LAN940 MUL118 MUL119 MUL121 MUL128 MUL129 MUL130 MUL413 MUL414 MUL416 MUL418 PAZ061 PAZ062 PAZ063 PAZ065 PAZ068 PAZ378 PAZ379 PAZ380 PAZ381 PAZ382 POR207 POR209 POR211 POR214 POR216 POR217 POR219 POR220 POR221 POR222 SOQ391 SOQ395 SOQ397 SOQ401 SOQ402 SOQ403 SOQ404 SOQ405 SOQ406 SOQ408 TOUR331 TOUR332 TOUR333 TOUR336 TOUR337 TOUR338 TOUR340 TOUR341 TOUR343 TOUR345;
do 
      bedtools coverage -a ../../00.input/bait.bed -b /dss/dsslegfs01/pr53da/pr53da-dss-0029/TranscriptomeAnalyses/1.Paleomix/${ind}.FullPaleoPipeline.bam -mean >> 01.output/${ind}_depth.bed

      bedtools coverage -a ../../00.input/bait.bed -b /dss/dsslegfs01/pr53da/pr53da-dss-0029/TranscriptomeAnalyses/1.Paleomix/${ind}.FullPaleoPipeline.bam -d >> 01.output/${ind}_depthperbp.bed
done

bedtools coverage -a ../../00.input/bait.bed -b /dss/dsslegfs01/pr53da/pr53da-dss-0029/TranscriptomeAnalyses/2.BedTools/MergingBams/MergedTranscriptomeBams.bam -mean >> 01.output/globdepth.bed

bedtools coverage -a ../../00.input/bait.bed -b /dss/dsslegfs01/pr53da/pr53da-dss-0029/TranscriptomeAnalyses/2.BedTools/MergingBams/MergedTranscriptomeBams.bam -d >> 01.output/globdepthbp.bed

ENDTIME=$(date +%s)
TIMESPEND=$(($ENDTIME - $STARTTIME))
((sec=TIMESPEND%60, TIMESPEND/=60, min=TIMESPEND%60, hrs=TIMESPEND/60))
timestamp=$(printf "%d:%02d:%02d" $hrs $min $sec)
echo "Job ended at $(date). Took $timestamp hours:minutes:seconds to complete."
```



### NGSadmix

To analyse the population structure, I used NGSadmix as implemented in ANGSD. NGSadmix takes genotype likelihood input data and assigns individuals to previously defined number of clusters K, based on maximising Hardy-Weinberg-equilibrium. 
I used the ``.beagle.gz`` file generated in [PCA](#pca) and a wrapper written by .... to specify the numbers of runs and Ks per job and submit them to slurm. I used conda `conda activate angsd933` because other versions of angsd don't work.  
To plot the admixture proportions per K, I used the following script in R:

```R
setwd("C:/Users/Jag/Documents/LMU/work/01.capture/02.ngsadmix")
colperpop <- c("GAB"="#001219", "HER"="#305762", "SOQ"="#447B8A", "TOU"="#79968E", "ARA"="#BEC6A5", "POR"="#FFF5BE", "MUL"="#F2D18A", "FOR"="#E4AD5B", "PAZ"="#CF7B4E", "LAN"="#C23436", "ESC"="#8D1216")


colors = c("#8D1216","#001219", #esc gab
           "#FFF5BE","#BEC6A5", #mul ara
           "#447B8A","#E4AD5B", #soq paz
           "#F2D18A","#305762", #tou her  
           "#CF7B4E","#79968E", #for tou,
           "grey") 

#gab her soq tou ara por mul for paz lan esc
mainlines = c(0,5,15,26,36,46,56,65,73,83,93,98)
dotlines = c(1:98)
# Colors you need to change manually,PITA
colorscrambles = list()
colorscrambles[[2]] = c(1,2)
colorscrambles[[3]] = c(2,1,3)
colorscrambles[[4]] = c(4,2,1,3)
colorscrambles[[5]] = c(3,2,4,5,1)
colorscrambles[[6]] = c(4,5,1,2,3,6)
colorscrambles[[7]] = c(4,3,5,6,2,1,7)
colorscrambles[[8]] = c(1,5,6,2,7,4,8,3)
colorscrambles[[9]] = c(1,3,2,4,9,8,7,6,5)
colorscrambles[[10]] = c(6,8,1,9,5,10,2,4,3,7)
colorscrambles[[11]] = c(7,4,11,10,8,6,5,2,1,9,3)

par(mar=c(0.75,0,0,0),oma=c(0,0,1,0))
layout(c(1:12)) #11 is max clusters

for (i in 2:11) { #11 is max clusters
  Q = t(as.matrix(read.table(paste0("beagle_bait_K",i,
                                    "_50reps.qopt"))))
  # The order of populations you want to show on the plot
  col.order <- c(24:28,29:38,79:88,56,89:98,1:10,69:78,49:55,57:58,16:23,59:68,39:48,11:15)
  Q <- Q[,col.order]
  barplot(Q,col=colors[colorscrambles[[i]]],border=colors[
    colorscrambles[[i]]],axes=F,space=c(rep(0,98)))
  sapply(mainlines,function(x){lines(c(x,x),c(0,1),lty=1,lwd=1.5)})
  sapply(dotlines,function(x){lines(c(x,x),c(0,1),lty=3,col="black")})
  mtext(paste0("K=",i),side=2,line=-2,font =2)
}
```

## F<sub>ST</sub> & IBD

### SAF
made saf files for each pop (MUL414 reassigned to SOQ414!) and ran on biohpc_gen partition using:
with `-minInd 4` change to 80% individuals
`-setMaxDepth` removed because pca same with and without

```bash
for pop in GAB HER SOQ TOU ARA POR MUL FOR PAZ LAN ESC; do
        # count lines in bamfile = number of individuals, extract only number)
        ind=$(wc -l /dss/dsslegfs01/pr53da/pr53da-dss-0029/capture/00.input/bamlist_${pop}.txt | sed -e 's/\s.*$//')
        echo "${ind}"
        # make 80% the min inds threshold 
        inds0=$(($ind * 8 ))
        inds=$(($inds0 / 10 | bc))
        echo "${inds}"
echo "#!/bin/bash
#SBATCH --clusters=biohpc_gen
#SBATCH --partition=biohpc_gen_normal
#SBATCH --mail-user=hagberg@lmu.de
#SBATCH --mail-type=ALL
#SBATCH --cpus-per-task=4
#SBATCH --mem-per-cpu=4000mb
#SBATCH --time=48:00:00
#SBATCH -J saf_${pop}
#SBATCH --out=saf_${pop}.out

angsd -b /dss/dsslegfs01/pr53da/pr53da-dss-0029/capture/00.input/bamlist_${pop}.txt -ref /dss/dsslegfs01/pr53da/pr53da-dss-0029/TranscriptomeAnalyses/TranscriptomeReference/grasshopperWolbRef.fasta -anc /dss/dsslegfs01/pr53da/pr53da-dss-0029/TranscriptomeAnalyses/TranscriptomeReference/grasshopperWolbRef.fasta -doSaf 1 -doCounts 1 -doMaf 1 -doMajorMinor 1 -GL 1 -r chr1 -sites /dss/dsslegfs01/pr53da/pr53da-dss-0029/capture/00.input/bait.sites -minInd ${inds} -minQ 20 -minMapQ 15 -only_proper_pairs 0 -remove_bads 1 -uniqueOnly 1 -C 50 -baq 1 -nThreads 6 -fold 0 -SNP_pval 1 -out saf_bait_${pop}" > 00.slurmscripts/saf_bait_${pop}.sh    

sbatch 00.slurmscripts/saf_bait_${pop}.sh

done
```

results in 0 sites => troubleshoot
Tries to solve:
1. module load angsd/0.933-gcc8 instead of conda (source to suggestion: [github question](https://github.com/ANGSD/angsd/issues/385)

> worked 
total: (60427/8897860 sites retained after filtering)

GAB | HER | SOQ | TOU | ARA | POR | MUL | FOR | PAZ | LAN | ESC
----|-----|-----|-----|-----|-----|-----|-----|-----|-----|-----|
60427 / 8897860 | 87517 / 9891295 | 96862 / 10295114 | 72717 / 10132607 | 68631 / 10162064 | 76725 / 10385723 | 78271 / 10070649 | 55741 / 9846189 | 55181 / 10038083 | 62874 / 10168991 | 48222 / 9022308

> how many sites should be retained though?
> this is a huge loss in sites, and will decrease further in the unlinked dataset. sad but nothing to do about it (quality>quantity)

I also made per individual SAF-files for downstream heterozygosity calculations using: 

```bash
#!/bin/bash
#SBATCH --clusters=biohpc_gen
#SBATCH --partition=biohpc_gen_normal
#SBATCH --mail-user=hagberg@lmu.de
#SBATCH --mail-type=ALL
#SBATCH --cpus-per-task=80
#SBATCH --mem-per-cpu=4000mb
#SBATCH --time=48:00:00
#SBATCH -J saf_ind
#SBATCH --out=saf_ind.out

module load angsd/0.933-gcc8

for ind in ARA270 ARA271 ARA273 ARA275 ARA276 ARA277 ARA279 ARA280 ARA281 ARA285 ESC011 ESC012 ESC013 ESC015 ESC344 FOR876 FOR877 FOR881 FOR884 FOR885 FOR887 FOR888 FOR890 GAB512 GAB513 GAB514 GAB515 GAB517 HER450 HER451 HER452 HER453 HER454 HER456 HER457 HER459 HER464 HER465 LAN927 LAN928 LAN929 LAN930 LAN931 LAN934 LAN935 LAN936 LAN938 LAN940 MUL118 MUL119 MUL121 MUL128 MUL129 MUL130 MUL413 MUL414 MUL416 MUL418 PAZ061 PAZ062 PAZ063 PAZ065 PAZ068 PAZ378 PAZ379 PAZ380 PAZ381 PAZ382 POR207 POR209 POR211 POR214 POR216 POR217 POR219 POR220 POR221 POR222 SOQ391 SOQ395 SOQ397 SOQ401 SOQ402 SOQ403 SOQ404 SOQ405 SOQ406 SOQ408 TOUR331 TOUR332 TOUR333 TOUR336 TOUR337 TOUR338 TOUR340 TOUR341 TOUR343 TOUR345;
do
        angsd -b /dss/dsslegfs01/pr53da/pr53da-dss-0029/capture/00.input/bam${ind}.txt -ref /dss/dsslegfs01/pr53da/pr53da-dss-0029/TranscriptomeAnalyses/TranscriptomeReference/grasshopperWolbRef.fasta -anc /dss/dsslegfs01/pr53da/pr53da-dss-0029/TranscriptomeAnalyses/TranscriptomeReference/grasshopperWolbRef.fasta -doSaf 1 
-doCounts 1 -doMaf 1 -doMajorMinor 1 -GL 1 -r chr1 -sites /dss/dsslegfs01/pr53da/pr53da-dss-0029/capture/00.input/bait.sites -minQ 20 -minMapQ 15 -only_proper_pairs 0 -remove_bads 1 -uniqueOnly 1 -C 50 -baq 1 -nThreads 160 -fold 0 -SNP_pval 1 -out 01.output/saf_bait_${ind}
done
```

### SFS
From the per population SAF files, I made one dimenional site frequency spectra (1dSFS) and two dimensional site frequency spectra (2dSFS) using realSFS as implemented in ANGSD using the following scripts. I also made per individual site frequency spectra (indsfs) for downstream analysis of heterozygosity.

1dSFS construction: 
```bash
for pop in GAB HER SOQ TOU ARA POR MUL FOR PAZ LAN ESC; do

echo "#!/bin/bash
#SBATCH --clusters=biohpc_gen
#SBATCH --partition=biohpc_gen_normal
#SBATCH --mail-user=hagberg@lmu.de
#SBATCH --mail-type=FAIL
#SBATCH --cpus-per-task=4
#SBATCH --mem-per-cpu=4000mb
#SBATCH --time=48:00:00
#SBATCH -J sfs_${pop}
#SBATCH --out=sfs_${pop}.out

module load angsd/0.933-gcc8

realSFS /dss/dsslegfs01/pr53da/pr53da-dss-0029/capture/saf/saf_bait_${pop}.saf.idx -P 4 > sfs_${pop}.sfs.em">00.slurmscripts/sfs_bait_${pop}.sh

sbatch 00.slurmscripts/sfs_bait_${pop}.sh
```

2dSFS construction: 
```bash
for pop1 in GAB HER SOQ TOU ARA POR MUL FOR PAZ LAN ESC; do
        for pop2 in GAB HER SOQ TOU ARA POR MUL FOR PAZ LAN ESC; do
                if [ ! -f "00.slurmscripts/2dsfs_${pop1}_${pop2}.sh" ] && [ ! -f "00.slurmscripts/2dsfs_${pop2}_${pop1}.sh" ] && [ ${pop1} != ${pop2} ] ; then    
echo "#!/bin/bash
#SBATCH --clusters=biohpc_gen
#SBATCH --partition=biohpc_gen_normal
#SBATCH --mail-user=hagberg@lmu.de
#SBATCH --mail-type=FAIL
#SBATCH --cpus-per-task=4
#SBATCH --mem-per-cpu=4000mb
#SBATCH --time=48:00:00
#SBATCH -J 2dsfs_${pop1}_${pop2}
#SBATCH --out=2dsfs_${pop1}_${pop2}.out

module load angsd/0.933-gcc8

realSFS /dss/dsslegfs01/pr53da/pr53da-dss-0029/capture/saf/01.output/saf_bait_${pop1}.saf.idx /dss/dsslegfs01/pr53da/pr53da-dss-0029/capture/saf/01.output/saf_bait_${pop2}.saf.idx -r chr1: -P 4 > 01.output/2dsfs_${pop1}_${pop2}.ml">00.slurmscripts/2dsfs_${pop1}_${pop2}.sh

sbatch 00.slurmscripts/2dsfs_${pop1}_${pop2}.sh

                fi
        done
done
```

Individual SFS construction: 
```bash
#!/bin/bash
#SBATCH --clusters=biohpc_gen
#SBATCH --partition=biohpc_gen_normal
#SBATCH --mail-user=hagberg@lmu.de
#SBATCH --mail-type=FAIL
#SBATCH --cpus-per-task=4
#SBATCH --mem-per-cpu=4000mb
#SBATCH --time=48:00:00
#SBATCH -J ind_sfs
#SBATCH --out=ind_sfs.out

module load angsd/0.933-gcc8

for ind in ARA270 ARA271 ARA273 ARA275 ARA276 ARA277 ARA279 ARA280 ARA281 ARA285 ESC011 ESC012 ESC013 ESC015 ESC344 FOR876 FOR877 FOR881 FOR884 FOR885 FOR887 FOR888 FOR890 GAB512 GAB513 GAB514 GAB515 GAB517 HER450 HER451 HER452 HER453 HER454 HER456 HER457 HER459 HER464 HER465 LAN927 LAN928 LAN929 LAN930 LAN931 LAN934 LAN935 LAN936 LAN938 LAN940 MUL118 MUL119 MUL121 MUL128 MUL129 MUL130 MUL413 MUL414 MUL416 MUL418 PAZ061 PAZ062 PAZ063 PAZ065 PAZ068 PAZ378 PAZ379 PAZ380 PAZ381 PAZ382 POR207 POR209 POR211 POR214 POR216 POR217 POR219 POR220 POR221 POR222 SOQ391 SOQ395 SOQ397 SOQ401 SOQ402 SOQ403 SOQ404 SOQ405 SOQ406 SOQ408 TOUR331 TOUR332 TOUR333 TOUR336 TOUR337 TOUR338 TOUR340 TOUR341 TOUR343 TOUR345;
do
        realSFS /dss/dsslegfs01/pr53da/pr53da-dss-0029/capture/saf/01.output/saf_bait_${ind}.saf.idx -r chr1 -P 160 > 01.outputs/${i}_HZ_sfs.ml

done
```

### F<sub>ST</sub>
I used the 2dSFS to calculate pairwise F<sub>ST</sub> between all populations. I wrote the following script to itirate the analysis over all combinations:

Global F<sub>ST</sub>

```bash
#!/bin/bash
#SBATCH --clusters=biohpc_gen
#SBATCH --partition=biohpc_gen_normal
#SBATCH --mail-user=hagberg@lmu.de
#SBATCH --mail-type=FAIL
#SBATCH --cpus-per-task=4
#SBATCH --mem-per-cpu=4000mb
#SBATCH --time=48:00:00
#SBATCH -J globfst_${pop1}_${pop2}
#SBATCH --out=globfst_${pop1}_${pop2}.out

module load angsd/0.933-gcc8

for pop1 in GAB HER SOQ TOU ARA POR MUL FOR PAZ LAN ESC; do
        for pop2 in GAB HER SOQ TOU ARA POR MUL FOR PAZ LAN ESC; do
                if [ ! -f "01.output/globfst_${pop1}_${pop2}*" ] && [ ! -f "01.output/globfst_${pop2}_${pop1}*" ] && [ ${pop1} != ${pop2} ] ; then

                        realSFS fst index /dss/dsslegfs01/pr53da/pr53da-dss-0029/capture/saf/01.output/saf_bait_${pop1}.saf.idx /dss/dsslegfs01/pr53da/pr53da-dss-0029/capture/saf/01.output/saf_bait_${pop2}.saf.idx -sfs ../sfs/01.output/2dsfs_${pop1}_${pop2}.ml -fstout 01.output/globfst_${pop1}_${pop2}
                else
                        echo "output containing ${pop1} and ${pop2} exists"

                fi
        done
done
```
Then i compiled the global weighted and unweighted FST values into one file.
```bash
#!/bin/bash
#SBATCH --clusters=biohpc_gen
#SBATCH --partition=biohpc_gen_normal
#SBATCH --mail-user=hagberg@lmu.de
#SBATCH --mail-type=FAIL
#SBATCH --cpus-per-task=4
#SBATCH --mem-per-cpu=4000mb
#SBATCH --time=48:00:00
#SBATCH -J comp_globfst
#SBATCH --out=comp_globfst.out

module load angsd/0.933-gcc8

for pop1 in GAB HER SOQ TOU ARA POR MUL FOR PAZ LAN ESC; do
        for pop2 in GAB HER SOQ TOU ARA POR MUL FOR PAZ LAN ESC; do
                if [ -f "01.output/globfst_${pop1}_${pop2}.fst.idx" ]
                then
                        x=$(realSFS fst stats /dss/dsslegfs01/pr53da/pr53da-dss-0029/capture/fst/01.output/globfst_${pop1}_${pop2}.fst.idx)
                        echo -e "${pop1}_${pop2}        ${x}" >> compiled_globfst_all.txt
                fi
        done
done
```
I plotted the global F<sub>ST</sub> in Isolation by Distance models in R using the following code:
```R
###Populationwise Fst

library(ade4)
library(ggplot2)
library(extrafont)
library(RColorBrewer)
library(reshape2)
library(tidyverse)
library(RColorBrewer)
library(viridis)
library(extrafont)

#Reading in files and creating the matrices
PopwiseFst <- read.csv(file = "baits_PopulationwiseFst.csv", header = T, row.names = 1, sep = ",", dec = ".")
PopwiseDist <- read.csv(file = "PopulationwiseDistance.csv", header = T, row.names = 1, sep = ",", dec = ".")

PopFstMatrix <- dist(PopwiseFst)
PopDistMatrix <- dist(PopwiseDist)

as.matrix(PopFstMatrix)
as.matrix(PopDistMatrix)

#Performing Mantel's test
mantel.rtest(PopFstMatrix, PopDistMatrix, nrepet = 1000)

#Scatterplot
PopFstXDistDF <- read.csv(file = "fstxdist.csv", header = F, row.names = 1, dec = ".", col.names = c("Pair","Fst","Distance.km"))


#Gray area is 0.95 confidence interval
ggplot(data = PopFstXDistDF, aes(x = Distance.km, y = Fst)) +
  geom_point() +
  geom_smooth(method='lm') +
  expand_limits(y=0) +
  scale_y_continuous(expand = c(0.01, 0.01), limits = c(0, 0.54)) +
  scale_x_continuous(expand = c(0.01, 0.01), limits = c(0.0, 21.1)) +
  ggtitle("Isolation by distance") +
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
```
Then i moved on to per site F<sub>ST</sub> calculations: 
```bash
realSFS fst print [pop1pop2].fst.idx > [pop1pop2].fst
```

The per gene F<sub>ST</sub> were calculated using the custom script `loopFst.pl`
per gene F<sub>ST</sub>:
```bash
perl loopFst.pl grasshopperRef.positions [pop1pop2].fst > genefst_[pop1pop2]
```

###  heterozygosity / fstat per region
Based on the per individual sfs, I first summed the output from each individual sfs into one `.ml` file, using: 


he did per site
loop through regions with bed file

## with GL repeat stats me (per region & per pop)
use correct filter (da fonseca)
pi theta, etc >100 sites per gene (2dsfs) angsd

### theta, dxy, taj. D 

### LD, rho
recombination

### allele freq. 

## Clinal shit (1 SNP/region)
make allele frequency file with relevant pops.. 

inbreeding cooefficint density ?
capturing introns

maxdepth
mindpth
restrict to new target coordinates of file

+mt wol

## References 
Schubert M, Ermini L, Sarkissian CD, Jónsson H, Ginolhac A, Schaefer R, Martin MD, Fernández R, Kircher M, McCue M, Willerslev E, and Orlando L. "Characterization of ancient and modern genomes by SNP detection and phylogenomic and metagenomic analysis using PALEOMIX". Nat Protoc. 2014 May;9(5):1056-82. doi: 10.1038/nprot.2014.063. Epub 2014 Apr 10. PubMed PMID: 24722405.