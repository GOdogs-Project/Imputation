#!/usr/bin/env bash

#! RUN : sbatch runImpute.sh 

#! sbatch directives begin here ###############################
#! How many whole nodes should be allocated?
#SBATCH --nodes=1
#! How many (MPI) tasks will there be in total? (<= nodes*32)
#! The skylake/skylake-himem nodes have 32 CPUs (cores) each.
#SBATCH --ntasks=8
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

PATH=$PATH:/rfs/project/rfs-x31eBTdMHgM/Software/local/bin/

LEN=`perl -lane 'if($F[1] eq "SN:'${CHR}'"){ $F[2] =~ s/LN://; print $F[2]; }' ${GENOME}.dict`
HAPS_FILE="${PROJECT}/shapeit/${PROJECT}.phased.impute.haplotypes"
LEGEND_FILE="${PROJECT}/shapeit/${PROJECT}.phased.impute.legend"
GWAS_HAPS_FILE="${GWAS_PLINK}/${GWAS_PLINK}.phased.haps"

impute2 -use_prephased_g -m ${MAPS}/chr${CHR}.cf3.1_map.txt -h ${HAPS_FILE} -l ${LEGEND_FILE} -known_haps_g ${GWAS_HAPS_FILE} -int 1 ${LEN} -allow_large_regions -Ne 200 -o ${GWAS_PLINK}.impute2 -phase
