for pop1 in  OUT GAB HER SOQ TOU ARA POR MUL TRO FOR PAZ LAN ESC STP AIN URD OTS ARI LEK MAR ALM VEN; do
#for pop1 in FOR; do
        for pop2 in OUT GAB HER SOQ TOU ARA POR MUL TRO FOR PAZ LAN ESC STP AIN URD OTS ARI LEK MAR ALM VEN; do
		if [ ! -f "00.slurmscripts/2dsfs_${pop1}_${pop2}.sh" ] && [ ! -f "00.slurmscripts/2dsfs_${pop2}_${pop1}.sh" ] && [ ${pop1} != ${pop2} ] ; then
echo "#!/bin/bash
#SBATCH --clusters=biohpc_gen
#SBATCH --partition=biohpc_gen_normal
#SBATCH --mail-user=hagberg@lmu.de
#SBATCH --mail-type=FAIL
#SBATCH --cpus-per-task=4
#SBATCH --mem-per-cpu=4000mb
#SBATCH --time=48:00:00
#SBATCH -J 2dsfs_${pop1}_${pop2}
#SBATCH --out=2dsfs_${pop1}_${pop2}.out

module load angsd/0.933-gcc8

realSFS /dss/dsslegfs01/pr53da/pr53da-dss-0029/capture/06.saf/01.output/saf_all_${pop1}.saf.idx /dss/dsslegfs01/pr53da/pr53da-dss-0029/capture/06.saf/01.output/saf_all_${pop2}.saf.idx -r chr1: -P 4 > 01.output/2dsfs_${pop1}_${pop2}.ml">00.slurmscripts/2dsfs_${pop1}_${pop2}.sh

sbatch 00.slurmscripts/2dsfs_${pop1}_${pop2}.sh

		fi
	done
done
