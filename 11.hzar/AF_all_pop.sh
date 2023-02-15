#!/bin/bash
#SBATCH --clusters=biohpc_gen
#SBATCH --partition=biohpc_gen_normal
#SBATCH --mail-user=hagberg@lmu.de
#SBATCH --mail-type=ALL
#SBATCH --cpus-per-task=4
#SBATCH --mem-per-cpu=4000mb
#SBATCH --time=48:00:00
#SBATCH -J af_all
#SBATCH --out=af_all.out

module load angsd/0.933-gcc8

for pop in GAB HER SOQ TOU ARA POR MUL TRO FOR PAZ LAN ESC; do

	angsd -b /dss/dsslegfs01/pr53da/pr53da-dss-0029/capture/00.input/bamlist_${pop}.txt -ref /dss/dsslegfs01/pr53da/pr53da-dss-0029/TranscriptomeAnalyses/TranscriptomeReference/grasshopperWolbRef.fasta -anc /dss/dsslegfs01/pr53da/pr53da-dss-0029/TranscriptomeAnalyses/TranscriptomeReference/grasshopperWolbRef.fasta -doMaf 1 -doMajorMinor 1 -GL 1 -r chr1 -sites /dss/dsslegfs01/pr53da/pr53da-dss-0029/capture/00.input/snps_richards.sites -nThreads 160 -only_proper_pairs 0 -SNP_pval 1 -out 01.output/try_hzar_${pop}

done
