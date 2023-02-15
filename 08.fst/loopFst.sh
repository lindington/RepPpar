#!/bin/bash
#SBATCH --clusters=biohpc_gen
#SBATCH --partition=biohpc_gen_normal
#SBATCH --mail-user=hagberg@lmu.de
#SBATCH --mail-type=FAIL
#SBATCH --cpus-per-task=4
#SBATCH --mem-per-cpu=4000mb
#SBATCH --time=48:00:00
#SBATCH -J gene_fst
#SBATCH --out=gene_fst.out

module load angsd/0.933-gcc8

for pop1 in GAB HER; do
        for pop2 in LAN ESC; do
		if [ -f "01.output/globfst_${pop1}_${pop2}.fst.idx" ]
		then
			realSFS fst print 01.output/globfst_${pop1}_${pop2}.fst.idx > 01.output/globfst_${pop1}_${pop2}.fst
            perl loopFst.pl ../00.input/bait.positions.csv 01.output/globfst_${pop1}_${pop2}.fst > 01.output/genefst_${pop1}_${pop2}
		fi
	done
done

for pop1 in STP AIN; do
        for pop2 in ALM VEN; do
		if [ -f "01.output/globfst_${pop1}_${pop2}.fst.idx" ]
		then
			realSFS fst print 01.output/globfst_${pop1}_${pop2}.fst.idx > 01.output/globfst_${pop1}_${pop2}.fst
            perl loopFst.pl ../00.input/bait.positions.csv 01.output/globfst_${pop1}_${pop2}.fst > 01.output/genefst_${pop1}_${pop2}
		fi
	done
done
