#!/bin/bash
#SBATCH --clusters=biohpc_gen
#SBATCH --partition=biohpc_gen_normal
#SBATCH --mail-user=hagberg@lmu.de
#SBATCH --mail-type=FAIL
#SBATCH --cpus-per-task=4
#SBATCH --mem-per-cpu=4000mb
#SBATCH --time=48:00:00
#SBATCH -J 2dsfs_MUL_TRO
#SBATCH --out=2dsfs_MUL_TRO.out

module load angsd/0.933-gcc8

realSFS /dss/dsslegfs01/pr53da/pr53da-dss-0029/capture/06.saf/01.output/saf_all_MUL.saf.idx /dss/dsslegfs01/pr53da/pr53da-dss-0029/capture/06.saf/01.output/saf_all_TRO.saf.idx -r chr1: -P 4 > 01.output/2dsfs_MUL_TRO.ml
