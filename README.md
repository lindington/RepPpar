# Repeatability of Introgression using Target/Capture sequencing in the Pyrenean <i>Pseudochorthippus parallelus</i> hybrid zone.

I contributed to the following components of this project: 
- Sequencing
  - Lab things
- Assembly
  - [Paleomix](/01.paleomix/)
  - [Depth](02.depth/)
- Population Structure: 
  - [PCA](/04.pca)
  - [NGSadmix](/05.ngsadmix/)
- Summary Stats
  - [F<sub>ST</sub>](/08.fst/)
  - [Heterozygosity](/09.hz/)
  - [π, Wattersons's θ, and Tajima's D](/10.pi_theta_taj/)
  - [d<sub>XY</sub>](/12.dxy/)
- Clines
  - [Geographic Clines](/11.hzar/)
  - [Genomic clines](/14.gghybrid)
- Demography
  - [dadi](13.dadi/)
> **NOTE**: individual MUL414 has to be reassigned to SOQ414 throughout the analysis! (shouldn't rename bamfiles)

## Assembly
#### Paleomix
The data used for this project came from two rounds of sequencing. The first was contained most individuals from the Portalet transect and an outgroup from Italy. The second round of sequencing contained all in individuals from the Basque transect and a few additional individuals from the Portalet transect.

I used the [BAM Pipeline](https://paleomix.readthedocs.io/en/stable/bam_pipeline/index.html) from Paleomix (Schubert et al., 2014) to create ``.bam`` files (compressed binary version of a SAM file, i.e. text-based format originally for storing biological sequences aligned to a reference sequence), on the basis of Enriques scripts from a previous project. See the relevant files in [01.paleomix](/01.paleomix/README.md).

## Depth

I calculated depth and coverage. See the relevant files in [02.depth](02.depth/).

#### Depth to decide which inds to exclude 

I used bedtools to calculate the global depth per baited region, global depth per base pair, depth per individual, and depth per individual per base pair: with the [full script](02.depth/bedtools/full_script.sh).

## Population Structure

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

To plot the PCA components, I ran the [pca plotting script](04.04.pca/bait_pca.R) locally on my laptop, because the R packages needed aren't installed on the server. 

#### NGSadmix

To analyse the population structure, I used NGSadmix as implemented in ANGSD. NGSadmix takes genotype likelihood input data and assigns individuals to previously defined number of clusters K, based on maximising Hardy-Weinberg-equilibrium. 

I used the ``.beagle.gz`` file generated in [PCA](#pca) and a wrapper written by .... to specify the numbers of runs and Ks per job and submit them to slurm. I used my a version of angsd i installed on conda (v. 0.933) because other versions of angsd don't work for this.
To plot the admixture proportions per K, I used the [ngsadmix plotting script](05.ngsadmix/ngsadmix.R) locally in R:

## Summary Stats

#### SAF
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

#### SFS
From the per population SAF files, I made one dimenional site frequency spectra ([1dSFS construction](07.sfs/1DSFS_all.sh) ) and two dimensional site frequency spectra ([2dSFS construction](07.sfs/2DSFS_all.sh)) using realSFS as implemented in ANGSD using the following scripts. I also made per individual site frequency spectra ([IndSFS construction](07.sfs/indSFS_all.sh) ) for downstream analysis of heterozygosity.

#### F<sub>ST</sub>
I used the 2dSFS to calculate pairwise F<sub>ST</sub> between all populations. I wrote the [Global F<sub>ST</sub> script](08.fst/FST_global_all.sh) to itirate the analysis over all combinations. Then i compiled the global weighted and unweighted FST values into one file using the [compiling script](08.fst/compile_FST_global_all.sh) I plotted the global F<sub>ST</sub> in Isolation by Distance models locally in R using the [global fst plotting script](08.fst/fst.R)

Then i moved on to per site F<sub>ST</sub> calculations: 
```bash
realSFS fst print [pop1pop2].fst.idx > [pop1pop2].fst
```

The per gene F<sub>ST</sub> were calculated using the custom script [loopFst.pl](08.fst/loopFst.pl) like this:
```bash
perl loopFst.pl grasshopperRef.positions [pop1pop2].fst > genefst_[pop1pop2]
```
I plotted the per gene fst using the [gene fest plotting script](08.fst/plot_gene_fst.R)

#### Heterozygosity
Based on the per individual sfs, I summed the output from each individual sfs into one `.ml` file, using the [summing script](09.hz/sum_indSFS.sh).

I plotted locally in R using the [hz plotting script](09.hz/hz.R).

#### Pi, Watterson's θ, Tajima's D (per region & per pop)
I calculated per site Watterson's θ, Tajima's D, and pi based on the 1DSFS using realSFS in angsd. I used the [site theta script](10.pi_theta_taj/site_thetas_all.sh) to calculate per site stats, then the [log theta script](10.pi_theta_taj/logthetas_all.sh) to extract the relevant output and then summed the sites into genes using the [gene thetas script](10.pi_theta_taj/gene_thetas.R) in R on the server by submitting it with the [gene theta bash script](10.pi_theta_taj/gene_thetas.sh). 

The output file contains many correlated summary stats. I extracted the relevant ones and plotted locally in R using my [sumstat plotting script](10.pi_theta_taj/plot_gene_theta.R).

#### d<sub>XY</sub>
I first calculated per site d<sub>XY</sub> using the the [site dxy script](12.dxy/calcDxy.R) written by Joshua Penalba, then summed into genes using [gene dxy script](12.dxy/gene_dxys.R) written by me and Alexander Hausmann. 

I plotted with the [d<sub>XY</sub> plotting script](12.plots.R).

## Clines
#### Geographic Clines

I made allele frequency files using many conversions to get from angsd compatible files to vcftools compatible files. The pipeline can be found in my [processing script](11.hzar/processing.sh).

I selected the SNPs using ... I ran hzar locally in R using my [geo clines script](11.hzar/geo_clines.R). I extracted relevant information using the [cline summary script](11.hzar/cline_summary.R) and plotted in python using [cline plotting script](11.hzar/clineplot.py) and calculated stats using [cline stats script](11.hzar/clinestats.py).

#### Genomic Clines
Richard Bailey and Maria conducted this Analysis, but i prepared the files.

I also made the Hybrid index geographic cline. I did this using the [hybrid index cline script](11.hzar/hicline.R).

## Demographic Analyses

#### dadi
Dörte Neumeister conducted this analysis but i prepared the files. 

## References 
Schubert M, Ermini L, Sarkissian CD, Jónsson H, Ginolhac A, Schaefer R, Martin MD, Fernández R, Kircher M, McCue M, Willerslev E, and Orlando L. "Characterization of ancient and modern genomes by SNP detection and phylogenomic and metagenomic analysis using PALEOMIX". Nat Protoc. 2014 May;9(5):1056-82. doi: 10.1038/nprot.2014.063. Epub 2014 Apr 10. PubMed PMID: 24722405.