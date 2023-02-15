#!/bin/bash
#SBATCH --clusters=biohpc_gen
#SBATCH --partition=biohpc_gen_normal
#SBATCH --mail-user=hagberg@lmu.de
#SBATCH --mail-type=ALL
#SBATCH --cpus-per-task=4
#SBATCH --mem-per-cpu=4000mb
#SBATCH --time=48:00:00
#SBATCH -J vcf_dadi
#SBATCH --out=vcf_dadi.out

module load angsd/0.933-gcc8

for l in linked random; do 
    for vs in parent all; do
        angsd -b /dss/dsslegfs01/pr53da/pr53da-dss-0029/capture/00.input/bamlist_dadi_${vs}.txt -ref /dss/dsslegfs01/pr53da/pr53da-dss-0029/TranscriptomeAnalyses/TranscriptomeReference/grasshopperWolbRef.fasta -anc /dss/dsslegfs01/pr53da/pr53da-dss-0029/TranscriptomeAnalyses/TranscriptomeReference/grasshopperWolbRef.fasta -doPlink 2 -doPost 1 -GL 1 -doMajorMinor 1 -doMaf 1 -doCounts 1 -doGeno -4 -r chr1 -sites /dss/dsslegfs01/pr53da/pr53da-dss-0029/capture/13.dadi/${l}snps_${vs}.sites -nThreads 160 -minQ 20 -minMapQ 15 -only_proper_pairs 0 -remove_bads 1 -uniqueOnly 1 -C 50 -baq 1 -SNP_pval 1 -out 01.output/vcf_dadi_${l}_${vs}
    done
done