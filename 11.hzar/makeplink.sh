#!/bin/bash
#SBATCH --clusters=biohpc_gen
#SBATCH --partition=biohpc_gen_normal
#SBATCH --mail-user=hagberg@lmu.de
#SBATCH --mail-type=ALL
#SBATCH --cpus-per-task=4
#SBATCH --mem-per-cpu=4000mb
#SBATCH --time=48:00:00
#SBATCH -J hzar_plink
#SBATCH --out=hzar_plink.out

module load angsd/0.933-gcc8

# count lines in bamfile = number of individuals, extract only number)
ind=$(wc -l /dss/dsslegfs01/pr53da/pr53da-dss-0029/capture/00.input/bamlist_bothtransects.txt | sed -e 's/\s.*$//')
echo "${ind}"
# make 80% the min inds threshold 
inds0=$(($ind * 8 ))
inds=$(($inds0 / 10 | bc))
echo "${inds}"

angsd -b /dss/dsslegfs01/pr53da/pr53da-dss-0029/capture/00.input/bamlist_bothtransects.txt -ref /dss/dsslegfs01/pr53da/pr53da-dss-0029/TranscriptomeAnalyses/TranscriptomeReference/grasshopperWolbRef.fasta -anc /dss/dsslegfs01/pr53da/pr53da-dss-0029/TranscriptomeAnalyses/TranscriptomeReference/grasshopperWolbRef.fasta -doPlink 2 -doPost 1 -GL 1 -doMajorMinor 1 -doMaf 1 -doCounts 1 -doGeno -4 -r chr1 -sites /dss/dsslegfs01/pr53da/pr53da-dss-0029/capture/13.dadi/linkedsnps_parent.sites -nThreads 160 -minInd ${inds} -fold 0 -SNP_pval 1e-2 -minQ 20 -minMapQ 15 -only_proper_pairs 0 -remove_bads 1 -uniqueOnly 1 -C 50 -baq 1 -out 01.output/hzar_plink_allinds