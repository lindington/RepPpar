#!/bin/bash
#SBATCH --clusters=biohpc_gen
#SBATCH --partition=biohpc_gen_normal
#SBATCH --mail-user=hagberg@lmu.de
#SBATCH --mail-type=FAIL
#SBATCH --cpus-per-task=4
#SBATCH --mem-per-cpu=4000mb
#SBATCH --time=48:00:00
#SBATCH -J comp_globfst
#SBATCH --out=comp_globfst.out

module load angsd/0.933-gcc8

for pop1 in OUT GAB HER SOQ TOU ARA POR MUL TRO FOR PAZ LAN ESC STP AIN URD OTS ARI LEK MAR ALM VEN; do
        for pop2 in OUT GAB HER SOQ TOU ARA POR MUL TRO FOR PAZ LAN ESC STP AIN URD OTS ARI LEK MAR ALM VEN; do
		if [ -f "01.output/globfst_${pop1}_${pop2}.fst.idx" ]
		then
			x=$(realSFS fst stats /dss/dsslegfs01/pr53da/pr53da-dss-0029/capture/08.fst/01.output/globfst_${pop1}_${pop2}.fst.idx)
			echo -e "${pop1}_${pop2}	${x}" >> compiled_globfst_all.txt
		fi
	done
done
