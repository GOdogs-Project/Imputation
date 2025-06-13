#!/usr/bin/env bash

#! RUN : sbatch addRefSNPs.sh 

#! sbatch directives begin here ###############################
#! How many whole nodes should be allocated?
#SBATCH --nodes=1
#! How many (MPI) tasks will there be in total? (<= nodes*32)
#! The skylake/skylake-himem nodes have 32 CPUs (cores) each.
#SBATCH --ntasks=1
#! How much wallclock time will be required?
#SBATCH --time=12:00:00
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

cd ${PROJECT}/addRefSNPs

VCF_ALL=${!PROJECT}
VCF=`dirname $VCF_ALL`
VCF+="/${PROJECT}-chr${CHR}.vcf.gz"

# remove split bed files
# rm -rf file*.bed

# concat all extra VCF files back to a single file - aka "GATHER"
ls extra.file*.gz > vcf-files.list
bcftools concat -f vcf-files.list | bgzip -c > ${PROJECT}.extra.vcf.gz

# call SNPs & remove them
bcftools call -Ov -mv ${PROJECT}.extra.vcf.gz > out.vcf

# create list of positions/SNPs to ignore in ref.uniq.snps
# remove positions with low read depth
grep -v '#' out.vcf | cut -f 3 > snps.list
bcftools query -e'FILTER="PASS"' -f'%ID\n' ${PROJECT}.extra.vcf.gz >> snps.list
bcftools view -H --exclude ID=@snps.list ${PROJECT}.extra.vcf.gz | cut -f 3 > ref.uniq.txt

rm -rf ${PROJECT}.extra.vcf.gz out.vcf snps.list

# call remainder of positions as ref/ref (0/0)
bcftools view -h ${VCF} | grep '#CHROM' > ${PROJECT}.extra.vcf
COUNT=`cut -f 10- ${PROJECT}.extra.vcf | tr "\t" "\n" | wc -l`
bcftools annotate -x FILTER ../../${PROJECT}.ref-panel.vcf.gz | bcftools view --include ID==@ref.uniq.txt -H -G | perl -lane 'print $_,"\tGT", "\t0/0"x'${COUNT}';' >> ${PROJECT}.extra.vcf

cp ../step1/ref-alleles .
bcftools annotate -x FILTER ../../${PROJECT}.ref-panel.vcf.gz | bcftools view --include ID==@ref.uniq.txt -H -G | cut -f 3,4 >> ref-alleles

# vcf2PLINK
eval "$PLINK --const-fid 0 --vcf ${PROJECT}.extra.vcf --out ${PROJECT}.extra"

# merge with .2 from step1
eval "$PLINK --bfile ../step1/${PROJECT}.2 --bmerge ${PROJECT}.extra --make-bed --out extra"
if [ -f "extra-merge.missnp" ]; then
  eval "$PLINK --bfile ${PROJECT}.extra --exclude extra-merge.missnp --make-bed --out temp"
  eval "$PLINK --bfile ../step1/${PROJECT}.2 --bmerge temp --make-bed --out ${PROJECT}.extra"
  # rm -rf temp* extra-merge*
fi

# filter SNPs and overwrite .3 from original run
# NB. do not use maf as all hom/ref AND hom/alt snps will be removed!
eval "$PLINK --bfile ${PROJECT}.extra --mind 0.1 --geno 0.03 --make-bed --out ${PROJECT}.3"
eval "$PLINK --bfile ${PROJECT}.3 --recode vcf-iid --a2-allele ref-alleles --real-ref-alleles --out ${PROJECT}"
bgzip ${PROJECT}.vcf; tabix -p vcf ${PROJECT}.vcf.gz

# export BED file of SNPs
awk -v OFS="\t" '{print $1, $4-1, $4, $2}' ${PROJECT}.3.bim > ${PROJECT}.snps.bed

rm -rf extra.* 
cp ${PROJECT}.snps.bed ../
cp ${PROJECT}.vcf.gz* ../