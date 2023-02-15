#!/bin/bash
#SBATCH -J 219_logtheta
#SBATCH --output=logThetas.out
#SBATCH --clusters=biohpc_gen
#SBATCH --partition=biohpc_gen_production
#SBATCH --cpus-per-task=40
#SBATCH --mem=600gb
#SBATCH --mail-type=ALL
#SBATCH --mail-user=hagberg@lmu.de
#SBATCH -t 4-00:00:00

module load angsd/0.933-gcc8

for pop in OUT GAB HER SOQ TOU ARA POR MUL TRO FOR PAZ LAN ESC STP AIN URD OTS ARI LEK MAR ALM VEN; do
    thetaStat print 01.output/site_theta_${pop}.thetas.idx > 01.output/logThetas_${pop}.chr1
done
