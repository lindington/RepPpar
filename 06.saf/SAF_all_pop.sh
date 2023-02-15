for pop in OUT GAB HER SOQ TOU ARA POR MUL TRO FOR PAZ LAN ESC STP AIN URD OTS ARI LEK MAR ALM VEN; do
#redo without low coverage ind in FOR (FOR881)
#for pop in FOR; do
	# count lines in bamfile = number of individuals, extract only number)
	ind=$(wc -l /dss/dsslegfs01/pr53da/pr53da-dss-0029/capture/00.input/bamlist_${pop}.txt | sed -e 's/\s.*$//')
	echo "${ind}"
	# make 80% the min inds threshold 
	inds0=$(($ind * 8 ))
	inds=$(($inds0 / 10 | bc))
	echo "${inds}"
echo "#!/bin/bash
#SBATCH --clusters=biohpc_gen
#SBATCH --partition=biohpc_gen_normal
#SBATCH --mail-user=hagberg@lmu.de
#SBATCH --mail-type=ALL
#SBATCH --cpus-per-task=4
#SBATCH --mem-per-cpu=4000mb
#SBATCH --time=48:00:00
#SBATCH -J saf_${pop}
#SBATCH --out=saf_${pop}.out

module load angsd/0.933-gcc8

angsd -b /dss/dsslegfs01/pr53da/pr53da-dss-0029/capture/00.input/bamlist_${pop}.txt -ref /dss/dsslegfs01/pr53da/pr53da-dss-0029/TranscriptomeAnalyses/TranscriptomeReference/grasshopperWolbRef.fasta -anc /dss/dsslegfs01/pr53da/pr53da-dss-0029/TranscriptomeAnalyses/TranscriptomeReference/grasshopperWolbRef.fasta -doSaf 1 -doCounts 1 -doMaf 1 -doMajorMinor 1 -GL 1 -r chr1 -sites /dss/dsslegfs01/pr53da/pr53da-dss-0029/capture/00.input/bait.sites -minInd ${inds} -minQ 20 -minMapQ 15 -only_proper_pairs 0 -remove_bads 1 -uniqueOnly 1 -C 50 -baq 1 -nThreads 6 -fold 0 -SNP_pval 1 -out 01.output/saf_all_${pop}" > 00.slurmscripts/saf_all_${pop}.sh

sbatch 00.slurmscripts/saf_all_${pop}.sh

done
