#! /bin/bash

#SBATCH -J paleomix
#SBATCH -o paleomix.out
#SBATCH --clusters=biohpc_gen
#SBATCH --partition=biohpc_gen_production
#SBATCH --mail-user=hagberg@lmu.de
#SBATCH --mail-type=ALL
#SBATCH --cpus-per-task=80
#SBATCH --mem-per-cpu=4000mb
#SBATCH --time=4-00:00:00

STARTTIME=$(date +%s)
echo " Job started at $(date)."

source /etc/profile.d/modules.sh

echo "
# Setup of User Spack
module use /dss/dsslegfs02/pn73se/pn73se-dss-0000/spack/modules/$LRZ_INSTRSET/linux*
module load user_spack
" >> $HOME/.bashrc

########################################################
# "$ conda activate paleomix" before running this file #
########################################################
module load /lrz/sys/spack/release/21.1.1/modules/x86_64/linux-sles15-x86_64/openjdk/11.0.2

source activate paleomix

paleomix bam run paleopipeline_new.yaml

ENDTIME=$(date +%s)
TIMESPEND=$(($ENDTIME - $STARTTIME))
((sec=TIMESPEND%60, TIMESPEND/=60, min=TIMESPEND%60, hrs=TIMESPEND/60))
timestamp=$(printf "%d:%02d:%02d" $hrs $min $sec)
echo "Job ended at $(date). Took $timestamp hours:minutes:seconds to complete."
