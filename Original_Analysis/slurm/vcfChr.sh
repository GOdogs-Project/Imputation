#!/usr/bin/env bash

#! RUN : sbatch vcfChr.sh <TAG> <CHR>

#! sbatch directives begin here ###############################
#! How many whole nodes should be allocated?
#SBATCH --nodes=1
#! How many (MPI) tasks will there be in total? (<= nodes*32)
#! The skylake/skylake-himem nodes have 32 CPUs (cores) each.
#SBATCH --ntasks=1
#! How much wallclock time will be required?
#SBATCH --time 2:00:00
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

module load bcftools-1.9-gcc-5.4.0-b2hdt5n        # bcftools
module load tabix-2013-12-16-gcc-5.4.0-xn3xiv7    # bgzip/tabix

TAG=$1
CHR=$2

source setup.config

VCF_ALL=${!TAG}
VCF=`dirname $VCF_ALL`
VCF+="/${TAG}-chr${CHR}.vcf.gz"

tabix -p vcf $VCF_ALL
if [ ! -e $VCF ]
then
  echo "Need to create $VCF";
  echo "chr$CHR $CHR" > $TAG.chrs
  bcftools view $VCF_ALL --regions $CHR,chr$CHR | bcftools annotate --rename-chrs $TAG.chrs | bcftools annotate --set-id +'%CHROM:%POS' | bgzip -c > $VCF
  tabix -p vcf $VCF
  rm -rf $TAG.chrs
fi
