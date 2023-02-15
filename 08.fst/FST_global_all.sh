#!/bin/bash
#SBATCH --clusters=biohpc_gen
#SBATCH --partition=biohpc_gen_normal
#SBATCH --mail-user=hagberg@lmu.de
#SBATCH --mail-type=FAIL
#SBATCH --cpus-per-task=4
#SBATCH --mem-per-cpu=4000mb
#SBATCH --time=48:00:00
#SBATCH -J globfst_all
#SBATCH --out=globfst_all.out

module load angsd/0.933-gcc8

#for pop1 in  OUT GAB HER SOQ TOU ARA POR MUL TRO FOR PAZ LAN ESC STP AIN URD OTS ARI LEK MAR ALM VEN; do
for pop1 in FOR; do
        for pop2 in OUT GAB HER SOQ TOU ARA POR MUL TRO FOR PAZ LAN ESC STP AIN URD OTS ARI LEK MAR ALM VEN; do
		if [ ! -f "01.output/globfst_${pop1}_${pop2}*" ] && [ ! -f "01.output/globfst_${pop2}_${pop1}*" ] && [ ${pop1} != ${pop2} ] ; then
			
			realSFS fst index /dss/dsslegfs01/pr53da/pr53da-dss-0029/capture/06.saf/01.output/saf_all_${pop1}.saf.idx /dss/dsslegfs01/pr53da/pr53da-dss-0029/capture/06.saf/01.output/saf_all_${pop2}.saf.idx -sfs ../07.sfs/01.output/2dsfs_${pop1}_${pop2}.ml -fstout 01.output/globfst_${pop1}_${pop2}
		else 
			echo "output containing ${pop1} and ${pop2} exists"

		fi
	done
done
