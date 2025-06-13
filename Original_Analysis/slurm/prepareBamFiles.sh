#!/usr/bin/env bash

#! RUN : sbatch prepareBamFiles.sh 

#! sbatch directives begin here ###############################
#! How many whole nodes should be allocated?
#SBATCH --nodes=1
#! How many (MPI) tasks will there be in total? (<= nodes*32)
#! The skylake/skylake-himem nodes have 32 CPUs (cores) each.
#SBATCH --ntasks=1
#! How much wallclock time will be required?
#SBATCH --time=0:30:00
#! What types of email messages do you wish to receive?
#SBATCH --mail-type=ALL
##SBATCH --mail-type=FAIL
#! Uncomment this to prevent the job from being requeued (e.g. if
#! interrupted by node failure or system downtime):
##SBATCH --no-requeue
#! For 6GB per CPU, set "-p skylake"; for 12GB per CPU, set "-p skylake-himem":
#SBATCH -p skylake

#SBATCH -o logs/job-%A_%a.out

module purge                               # Removes all modules still loaded
module load rhel7/default-peta4            # REQUIRED - loads the basic environment

module load samtools/1.9                          # samtools

PROJECT=$1
CHR=$2
source setup.config

FILE=$(ls ${BAM_FILES}/*.bam | sed -n ${SLURM_ARRAY_TASK_ID}p)
F_NAME=`basename ${FILE}`
samtools view -b ${FILE} ${CHR} > ${BAM_FILES}/chr${CHR}/chr${CHR}-${F_NAME}
samtools index ${BAM_FILES}/chr${CHR}/chr${CHR}-${F_NAME}
