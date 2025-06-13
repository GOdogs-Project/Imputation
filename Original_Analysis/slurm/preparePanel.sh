#!/usr/bin/env bash

#! RUN : sbatch preparePanel.sh 

#! sbatch directives begin here ###############################
#! How many whole nodes should be allocated?
#SBATCH --nodes=1
#! How many (MPI) tasks will there be in total? (<= nodes*32)
#! The skylake/skylake-himem nodes have 32 CPUs (cores) each.
#SBATCH --ntasks=1
#! How much wallclock time will be required?
#SBATCH --time=00:30:00
#! What types of email messages do you wish to receive?
#SBATCH --mail-type=ALL
##SBATCH --mail-type=FAIL
#! Uncomment this to prevent the job from being requeued (e.g. if
#! interrupted by node failure or system downtime):
##SBATCH --no-requeue
#! For 6GB per CPU, set "-p skylake"; for 12GB per CPU, set "-p skylake-himem":
#SBATCH -p skylake

#SBATCH -o logs/job-%j.out

module purge                               # Removes all modules still loaded
module load rhel7/default-peta4            # REQUIRED - loads the basic environment

module load bcftools-1.9-gcc-5.4.0-b2hdt5n        # bcftools
module load tabix-2013-12-16-gcc-5.4.0-xn3xiv7    # bgzip/tabix
module load plink-1.9-gcc-5.4.0-sm3ojoi           # plink/plinkdog

PROJECT=$1
CHR=$2
source setup.config

cd ${PROJECT}
mkdir shapeit; cd shapeit

bcftools merge -m id --force-samples --missing-to-ref ../../${PROJECT}.ref-panel.vcf.gz ../${PROJECT}.vcf.gz > merged.vcf
bcftools query -l merged.vcf| grep ':' > samples.dups
bcftools view -Ov -S ^samples.dups merged.vcf | bcftools view -m2 -M2 -v snps | bgzip -c > ${PROJECT}.panel.vcf.gz
tabix -p vcf ${PROJECT}.panel.vcf.gz

# mv ${PROJECT}.panel.vcf.gz* ../
rm -rf merged.vcf samples.dups

eval "$PLINK --const-fid 0 --vcf ${PROJECT}.panel.vcf.gz --make-bed"
perl -lane 'if(($F[4] eq "A" && $F[5] eq "T") || ($F[4] eq "T" && $F[5] eq "A") || ($F[4] eq "C" && $F[5] eq "G") || ($F[4] eq "G" && $F[5] eq "C")){ print $F[1] }' plink.bim > ambiguous.snps
eval "$PLINK --const-fid 0 --vcf ${PROJECT}.panel.vcf.gz --exclude ambiguous.snps --freq --missing --make-bed"
eval "$PLINK --bfile plink --maf 0.01 --mind 0.1 --geno 0.03 --make-bed --out ${PROJECT}.ref-panel"



