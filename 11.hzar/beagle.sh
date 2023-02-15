#! /bin/bash

#SBATCH -J 219_beagle
#SBATCH --output=219_beagle.out
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

angsd -b /dss/dsslegfs01/pr53da/pr53da-dss-0029/capture/00.input/bamlist_out.txt -ref /dss/dsslegfs01/pr53da/pr53da-dss-0029/TranscriptomeAnalyses/TranscriptomeReference/grasshopperWolbRef.fasta -doMajorMinor 1 -GL 1 -doGlf 2 -SNP_pval 1e-2 -doMaf 1 -nThreads 4 -r chr1: -sites /dss/dsslegfs01/pr53da/pr53da-dss-0029/capture/00.input/bait.sites -baq 1 -remove_bads 1 -uniqueOnly 1 -C 50 -minMapQ 15 -only_proper_pairs 0 -minQ 20 -doCounts 1 -doPost 2 -doGeno 32 -minInd 175 -setMaxDepth 14472 -out 01.output/219_out_maxdepth

ENDTIME=$(date +%s)
TIMESPEND=$(($ENDTIME - $STARTTIME))
((sec=TIMESPEND%60, TIMESPEND/=60, min=TIMESPEND%60, hrs=TIMESPEND/60))
timestamp=$(printf "%d:%02d:%02d" $hrs $min $sec)
echo "Job ended at $(date). Took $timestamp hours:minutes:seconds to complete."