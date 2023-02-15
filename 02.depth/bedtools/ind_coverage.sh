#! /bin/bash

#SBATCH --clusters=biohpc_gen
#SBATCH --partition=biohpc_gen_normal
#SBATCH --mail-user=hagberg@lmu.de
#SBATCH --mail-type=FAIL
#SBATCH --cpus-per-task=10
#SBATCH --mem-per-cpu=2000mb
#SBATCH --time=2-00:00:00
#SBATCH -J indcov
#SBATCH -o ind_cov.out

STARTTIME=$(date +%s)
echo " ###### Job name: $SLURM_JOB_NAME, job ID: $SLURM_JOB_ID ######"     
echo " Job started at $(date)."

module load bedtools2

for ind in RP672 RP674 RP675 RP676 RP677 RP679 RP680 RP682 RP684 RP687 GAB512 GAB513 GAB514 GAB515 GAB517 HER450 HER451 HER452 HER453 HER454 HER456 HER457 HER459 HER464 HER465 SOQ391 SOQ395 SOQ397 SOQ401 SOQ402 SOQ403 SOQ404 SOQ405 SOQ406 SOQ408 TOUR331 TOUR332 TOUR333 TOUR336 TOUR337 TOUR338 TOUR340 TOUR341 TOUR343 TOUR345 ARA270 ARA271 ARA273 ARA275 ARA276 ARA277 ARA279 ARA280 ARA281 ARA285 POR207 POR209 POR211 POR214 POR216 POR217 POR219 POR220 POR221 POR222 POR240 POR241 POR242 MUL118 MUL119 MUL121 MUL128 MUL129 MUL130 MUL413 MUL414 MUL416 MUL418 MUL123 MUL124 MUL125 MUL126 TRO170 TRO171 TRO172 TRO173 TRO174 TRO177 TRO178 TRO192 TRO193 TRO194 TRO196 TRO198 TRO200 FOR876 FOR877 FOR881 FOR884 FOR885 FOR887 FOR888 FOR890 FOR883 FOR896 FOR897 PAZ061 PAZ062 PAZ063 PAZ065 PAZ068 PAZ378 PAZ379 PAZ380 PAZ381 PAZ382 LAN927 LAN928 LAN929 LAN930 LAN931 LAN934 LAN935 LAN936 LAN938 LAN940 ESC011 ESC012 ESC013 ESC015 ESC344 PE261 PE262 PE268 PE276 PE278 PE279 PE280 PE281 PE282 PE283 AI149 AI151 AI152 AI170 AI173 AI178 AI179 AI180 AI183 AI185 UR18 UR26 UR27 UR28 UR3 UR31 UR4 UR5 UR8 UR9 OTS838 OTS839 OTS844 OTS845 OTS847 OTS850 OTS851 OTS852 OTS854 OTS855 ARI1001 ARI1002 ARI1004 ARI989 ARI992 ARI993 ARI994 ARI996 ARI997 ARI998 LEK811 LEK813 LEK814 LEK815 LEK816 LEK826 LEK827 LEK828 LEK829 MAR785 MAR789 MAR790 MAR791 MAR792 MAR799 MAR800 MAR801 MAR806 MAR807 ALM729 ALM730 ALM732 ALM734 ALM735 ALM736 ALM737 ALM738 ALM740 ALM741 VEN699 VEN700 VEN701 VEN702 VEN703 VEN704 VEN707 VEN708 VEN709 VEN710; do
	if [ -f "/dss/dsslegfs01/pr53da/pr53da-dss-0029/capture/01.paleomix/${ind}.FullPaleoPipeline.bam" ]; then
		bedtools coverage -a /dss/dsslegfs01/pr53da/pr53da-dss-0029/capture/00.input/bait.bed -b /dss/dsslegfs01/pr53da/pr53da-dss-0029/capture/01.paleomix/${ind}.FullPaleoPipeline.bam -mean > 01.output/ind_gcov_${ind}.bed
		#bedtools coverage -a /dss/dsslegfs01/pr53da/pr53da-dss-0029/capture/00.input/bait.bed -b /dss/dsslegfs01/pr53da/pr53da-dss-0029/capture/01.paleomix/${ind}.FullPaleoPipeline.bam -d > 01.output/ind_bpcov_${ind}.bed
	elif [ -f "/dss/dsslegfs01/pr53da/pr53da-dss-0029/TranscriptomeAnalyses/1.Paleomix/${ind}.FullPaleoPipeline.bam" ]; then
		bedtools coverage -a /dss/dsslegfs01/pr53da/pr53da-dss-0029/capture/00.input/bait.bed -b /dss/dsslegfs01/pr53da/pr53da-dss-0029/TranscriptomeAnalyses/1.Paleomix/${ind}.FullPaleoPipeline.bam -mean > 01.output/ind_gcov_${ind}.bed
		#bedtools coverage -a /dss/dsslegfs01/pr53da/pr53da-dss-0029/capture/00.input/bait.bed -b /dss/dsslegfs01/pr53da/pr53da-dss-0029/TranscriptomeAnalyses/1.Paleomix/${ind}.FullPaleoPipeline.bam -d > 01.output/ind_bpcov_${ind}.bed
	fi
done

ENDTIME=$(date +%s)
TIMESPEND=$(($ENDTIME - $STARTTIME))
((sec=TIMESPEND%60, TIMESPEND/=60, min=TIMESPEND%60, hrs=TIMESPEND/60))
timestamp=$(printf "%d:%02d:%02d" $hrs $min $sec)
echo "Job ended at $(date). Took $timestamp hours:minutes:seconds to complete."
