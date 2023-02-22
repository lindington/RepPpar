# Repeatability of Introgression using Target/Capture sequencing in the Pyrenean <i>Pseudochorthippus parallelus</i> hybrid zone.

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

> **NOTE**: individual MUL414 has to be reassigned to SOQ414 throughout the analysis! (shouldn't rename bamfiles)

## Assembly
#### Paleomix
The data used for this project came from two rounds of sequencing. The first was contained most individuals from the Portalet transect and an outgroup from Italy. The second round of sequencing contained all in individuals from the Basque transect and a few additional individuals from the Portalet transect.

I used the [BAM Pipeline](https://paleomix.readthedocs.io/en/stable/bam_pipeline/index.html) from Paleomix (Schubert et al., 2014) to create ``.bam`` files (compressed binary version of a SAM file, i.e. text-based format originally for storing biological sequences aligned to a reference sequence), on the basis of Enriques scripts from a previous project. See the relevant files in [01.paleomix](/01.paleomix/README.md).

## Depth

I calculated depth and coverage. See the relevant files in [02.depth](02.depth/).

### Depth to decide which inds to exclude 

I used bedtools to calculate the global depth per baited region, global depth per base pair, depth per individual, and depth per individual per base pair: with the [full script](02.depth/bedtools/full_script.sh).

## PCA and NGSadmix

For the population structure analyses, I used the full dataset including a geographical "outgroup" (bamlist_out.txt). To do any analysis in angsd, I needed a bamlist (`.txt` files containing paths to relevant bamfiles), which i made using `readlink -f ../01.paleomix/*bam > ../00.input/bamlist122.txt`.  

I started by making a `.geno` and a `.beagle` file with the same filters. The ``.beagle`` file will be used for [NGSadmix](#ngsadmix). I specified to retain the baited regions only used the `-sites` flag and restricted analyses to chr1 using the `-r` flag.

I ran [this script](03.beagle/beagle.sh) on the module installation of angsd (v. 0.933-gcc8) to make `.beagle` and `.geno` files.

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

To plot the PCA components, I ran the [plotting script](04.04.pca/bait_pca.R) locally on my laptop, because the R packages needed aren't installed on the server. 

### NGSadmix

To analyse the population structure, I used NGSadmix as implemented in ANGSD. NGSadmix takes genotype likelihood input data and assigns individuals to previously defined number of clusters K, based on maximising Hardy-Weinberg-equilibrium. 

I used the ``.beagle.gz`` file generated in [PCA](#pca) and a wrapper written by .... to specify the numbers of runs and Ks per job and submit them to slurm. I used my a version of angsd i installed on conda (v. 0.933) because other versions of angsd don't work for this.
To plot the admixture proportions per K, I used the [plotting script](05.ngsadmix/ngsadmix.R) locally in R:

## F<sub>ST</sub> & IBD

### SAF
I made ``.saf`` files for each pop (MUL414 reassigned to SOQ414!) by running the [per-pop saf script](06.saf/SAF_all_pop.sh) on biohpc_gen partition with: 
- `-minInd 4` changed to 80% individuals
- `-setMaxDepth` removed because pca not affected by it (same with and without)

> **NOTE**: conda version of angsd did not weork (results in 0 sites). module load angsd/0.933-gcc8 instead of conda (source to suggestion: [github question](https://github.com/ANGSD/angsd/issues/385))

total: (60427/8897860 sites retained after filtering)

GAB | HER | SOQ | TOU | ARA | POR | MUL | FOR | PAZ | LAN | ESC
----|-----|-----|-----|-----|-----|-----|-----|-----|-----|-----|
60427 / 8897860 | 87517 / 9891295 | 96862 / 10295114 | 72717 / 10132607 | 68631 / 10162064 | 76725 / 10385723 | 78271 / 10070649 | 55741 / 9846189 | 55181 / 10038083 | 62874 / 10168991 | 48222 / 9022308

> how many sites should be retained?
> this is a huge loss in sites, and will decrease further in the unlinked dataset. sad but nothing to do about it (quality>quantity)

I also made per individual SAF-files for downstream heterozygosity calculations using the [per-ind-saf script](06.saf/SAF_all_ind.sh).

### SFS
From the per population SAF files, I made one dimenional site frequency spectra ([1dSFS construction](07.sfs/1DSFS_all.sh) ) and two dimensional site frequency spectra ([2dSFS construction](07.sfs/2DSFS_all.sh)) using realSFS as implemented in ANGSD using the following scripts. I also made per individual site frequency spectra ([IndSFS construction](07.sfs/indSFS_all.sh) ) for downstream analysis of heterozygosity.

### F<sub>ST</sub>
I used the 2dSFS to calculate pairwise F<sub>ST</sub> between all populations. I wrote the [Global F<sub>ST</sub> script](08.fst/FST_global_all.sh) to itirate the analysis over all combinations. 

Then i compiled the global weighted and unweighted FST values into one file using the [compiling script](08.fst/compile_FST_global_all.sh)

I plotted the global F<sub>ST</sub> in Isolation by Distance models locally in R using the [plotting script](08.fst/fst.R)

Then i moved on to per site F<sub>ST</sub> calculations: 
```bash
realSFS fst print [pop1pop2].fst.idx > [pop1pop2].fst
```

The per gene F<sub>ST</sub> were calculated using the custom script [loopFst.pl](08.fst/loopFst.pl) like this:
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