#for pop in OUT GAB HER SOQ TOU ARA POR MUL TRO FOR PAZ LAN ESC STP AIN URD OTS ARI LEK MAR ALM VEN; do
for pop in FOR; do
echo "#!/bin/bash
#SBATCH --clusters=biohpc_gen
#SBATCH --partition=biohpc_gen_normal
#SBATCH --mail-user=hagberg@lmu.de
#SBATCH --mail-type=FAIL
#SBATCH --cpus-per-task=4
#SBATCH --mem-per-cpu=4000mb
#SBATCH --time=48:00:00
#SBATCH -J sfs_${pop}
#SBATCH --out=sfs_${pop}.out

module load angsd/0.933-gcc8

realSFS /dss/dsslegfs01/pr53da/pr53da-dss-0029/capture/06.saf/01.output/saf_all_${pop}.saf.idx -P 4 > 01.output/sfs_${pop}.sfs.em">00.slurmscripts/sfs_all_${pop}.sh

sbatch 00.slurmscripts/sfs_all_${pop}.sh

done
