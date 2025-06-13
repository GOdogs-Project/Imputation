#!/usr/bin/env bash

#! RUN : sbatch shapeitPanel.sh 

#! sbatch directives begin here ###############################
#! How many whole nodes should be allocated?
#SBATCH --nodes=1
#! How many (MPI) tasks will there be in total? (<= nodes*32)
#! The skylake/skylake-himem nodes have 32 CPUs (cores) each.
#SBATCH --ntasks=20
#! How much wallclock time will be required?
#SBATCH --time 12:00:00
#! What types of email messages do you wish to receive?
#SBATCH --mail-type=ALL
#! Uncomment this to prevent the job from being requeued (e.g. if
#! interrupted by node failure or system downtime):
##SBATCH --no-requeue

#! For 6GB per CPU, set "-p skylake"; for 12GB per CPU, set "-p skylake-himem":
#SBATCH -p skylake

#SBATCH -o logs/job-%j.out

module purge                               # Removes all modules still loaded
module load rhel7/default-peta4            # REQUIRED - loads the basic environment

PROJECT=$1
CHR=$2
source setup.config

mkdir ${PROJECT}/shapeit
cd ${PROJECT}/shapeit

PATH=$PATH:/rfs/project/rfs-x31eBTdMHgM/Software/local/bin/

shapeit -M ${MAPS}/chr${CHR}.cf3.1_map.txt -B ${PROJECT}.ref-panel -O ${PROJECT}.phased -T 20 --window 2 --effective-size 200 --force
shapeit -convert --input-haps ${PROJECT}.phased --output-ref ${PROJECT}.phased.impute
