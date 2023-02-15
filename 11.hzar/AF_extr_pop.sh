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

for pop in PAR ERY; do

	# count lines in bamfile = number of individuals, extract only number)
	ind=$(wc -l /dss/dsslegfs01/pr53da/pr53da-dss-0029/capture/00.input/bamlist_cline${pop}.txt | sed -e 's/\s.*$//')
	echo "${ind}"
	# make 60% the min inds threshold 
	inds0=$(($ind * 6 ))
	inds=$(($inds0 / 10 | bc))
	echo "${inds}"

	angsd -b /dss/dsslegfs01/pr53da/pr53da-dss-0029/capture/00.input/bamlist_cline${pop}.txt -ref /dss/dsslegfs01/pr53da/pr53da-dss-0029/TranscriptomeAnalyses/TranscriptomeReference/grasshopperWolbRef.fasta -anc /dss/dsslegfs01/pr53da/pr53da-dss-0029/TranscriptomeAnalyses/TranscriptomeReference/grasshopperWolbRef.fasta -doMaf 1 -doMajorMinor 1 -GL 1 -r chr1 -sites /dss/dsslegfs01/pr53da/pr53da-dss-0029/capture/00.input/bait.sites -minInd ${inds} -minQ 20 -minMapQ 15 -only_proper_pairs 0 -remove_bads 1 -uniqueOnly 1 -C 50 -baq 1 -nThreads 6 -fold 0 -dobcf 1 -dopost 1 -SNP_pval 1e-2 --ignore-RG 0 -dogeno 1 -docounts 1 -out 01.output/af_cline_${pop}

done
