#!/usr/bin/env bash

#! RUN : bash prepGwasData_wrapper.sh <CHR>
#! Eg. : bash prepGwasData_wrapper.sh 13

CHR=$1
SCRIPTS=`dirname $0`

[[ -z "$CHR" ]] && { echo "ERROR: No CHR provided for this run"; exit 1; }

mkdir -p chr${CHR}/logs; cd chr${CHR};

cp ../setup.config .
source setup.config

mkdir $GWAS_PLINK; cd $GWAS_PLINK

jid1=$(sbatch -A ${ACCOUNT} -J ${GWAS_PLINK}.prep ${SCRIPTS}/slurm/prepareGWAS.sh ${CHR})
sbatch -A ${ACCOUNT} -J ${GWAS_PLINK}.shapeit --dependency=afterok:${jid1##* } ${SCRIPTS}/slurm/shapeitGWAS.sh ${CHR}