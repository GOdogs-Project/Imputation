#!/usr/bin/env bash

#! RUN : sbatch shapeitGWAS.sh 

#! sbatch directives begin here ###############################
#! How many whole nodes should be allocated?
#SBATCH --nodes=1
#! How many (MPI) tasks will there be in total? (<= nodes*32)
#! The skylake/skylake-himem nodes have 32 CPUs (cores) each.
#SBATCH --ntasks=8
#! How much wallclock time will be required?
#SBATCH --time 00:30:00
#! What types of email messages do you wish to receive?
#SBATCH --mail-type=ALL
#! Uncomment this to prevent the job from being requeued (e.g. if
#! interrupted by node failure or system downtime):
##SBATCH --no-requeue

#! For 6GB per CPU, set "-p skylake"; for 12GB per CPU, set "-p skylake-himem":
#SBATCH -p skylake

#SBATCH -o ../logs/job-%j.out

module purge                               # Removes all modules still loaded
module load rhel7/default-peta4            # REQUIRED - loads the basic environment

CHR=$1
source ../setup.config

PATH=$PATH:/rfs/project/rfs-x31eBTdMHgM/Software/local/bin/

shapeit -M ${MAPS}/chr${CHR}.cf3.1_map.txt -B ${GWAS_PLINK} -O ${GWAS_PLINK}.phased -T 8 --window 2 --effective-size 200 --force
shapeit -convert --input-haps ${GWAS_PLINK}.phased --output-ref ${GWAS_PLINK}.phased.impute
