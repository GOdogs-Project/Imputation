#!/usr/bin/env bash

#! RUN : bash imputaion_wrapper.sh <CHR> <PROJECT>
#! Eg. : bash imputaion_wrapper.sh 13 go_lab

CHR=$1
PROJECT=$2
SCRIPTS=`dirname $0`

[[ -z "$CHR" ]] && { echo "ERROR: No CHR provided for this run"; exit 1; }
[[ -z "$PROJECT" ]] && { echo "ERROR: No PROJECT provided for this run"; exit 1; }
[[ ! -e chr${CHR} ]] && { echo "No directory exists for chr ${CHR}"; exit 1; }

cd chr${CHR};
source setup.config


# Check that the PROJECT haplotypes and legend files exist
if [[ ! -e ${PROJECT}/shapeit/${PROJECT}.phased.impute.haplotypes && -e ${PROJECT}/shapeit/${PROJECT}.phased.impute.legend ]]
then 
  echo "ERROR - Unable to find haplotypes/legend files for ${PROJECT} - have you run generatePanel_wrapper.sh first?";
  exit 1;
fi

# Check that the known haplotypes from teh GWAS data exist
if [ ! -e ${GWAS_PLINK}/${GWAS_PLINK}.phased.haps ]
then 
  echo "ERROR - Unable to find GWAS known haplotypes for ${GWAS_PLINK} - have you run prepGwasData_wrapper.sh first?";
  exit 1;
fi

# Files exist so we are ready to impute!
sbatch -A ${ACCOUNT} -J ${GWAS_PLINK}.impute ${SCRIPTS}/slurm/runImpute.sh ${PROJECT} ${CHR}
