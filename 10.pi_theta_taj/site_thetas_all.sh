#!/bin/bash
#SBATCH -J site_theta
#SBATCH --output=site_theta.out
#SBATCH --clusters=biohpc_gen
#SBATCH --partition=biohpc_gen_production
#SBATCH --cpus-per-task=40
#SBATCH --mem=600gb
#SBATCH --mail-type=ALL
#SBATCH --mail-user=hagberg@lmu.de
#SBATCH -t 4-00:00:00

module load angsd/0.933-gcc8

for pop in OUT GAB HER SOQ TOU ARA POR MUL TRO FOR PAZ LAN ESC STP AIN URD OTS ARI LEK MAR ALM VEN; do

realSFS saf2theta ../06.saf/01.output/saf_all_${pop}.saf.idx -sfs ../07.sfs/01.output/sfs_${pop}.sfs.em -outname 01.output/site_theta_${pop}

done