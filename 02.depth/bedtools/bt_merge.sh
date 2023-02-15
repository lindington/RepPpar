#! /bin/bash

#SBATCH --clusters=biohpc_gen
#SBATCH --partition=biohpc_gen_normal
#SBATCH --mail-user=hagberg@lmu.de
#SBATCH --mail-type=FAIL
#SBATCH --cpus-per-task=10
#SBATCH --mem-per-cpu=2000mb
#SBATCH --time=2-00:00:00
#SBATCH -J merge
#SBATCH -o merged.out

STARTTIME=$(date +%s)

echo " ###### Job name: $SLURM_JOB_NAME, job ID: $SLURM_JOB_ID ######"     
echo " Job started at $(date)."

bedtools merge -i /dss/dsslegfs01/pr53da/pr53da-dss-0029/capture/02.depth/bedtools/compilebeds_sorted.bed > /dss/dsslegfs01/pr53da/pr53da-dss-0029/capture/02.depth/bedtools/compilebeds_sorted_merged.bed

ENDTIME=$(date +%s)
TIMESPEND=$(($ENDTIME - $STARTTIME))
((sec=TIMESPEND%60, TIMESPEND/=60, min=TIMESPEND%60, hrs=TIMESPEND/60))    
timestamp=$(printf "%d:%02d:%02d" $hrs $min $sec)
echo "Job ended at $(date). Took $timestamp hours:minutes:seconds to complete."
