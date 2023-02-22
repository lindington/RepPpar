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