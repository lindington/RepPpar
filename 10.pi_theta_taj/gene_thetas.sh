#!/bin/bash
#SBATCH -J gene_theta
#SBATCH --output=gene_theta.out
#SBATCH --clusters=biohpc_gen
#SBATCH --partition=biohpc_gen_production
#SBATCH --cpus-per-task=40
#SBATCH --mem=600gb
#SBATCH --mail-type=ALL
#SBATCH --mail-user=hagberg@lmu.de
#SBATCH -t 4-00:00:00

module load gcc/6.5.0
module load r/4.0.5-gcc11-mkl
# no clue why gcc is needed suddenly, but it didnt work without it.

Rscript gene_thetas.R