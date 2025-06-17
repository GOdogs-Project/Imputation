# Imputation Analysis for Eleanor Raffan Lab
Analysis carried out by Anton Enright (aje39@cam.ac.uk) in consultation with Jade Scardham and Eleanor Raffan (PDN). The workflow used here is adapted from the GODogs Imputation workflow elsewhere on github.
This is part of the [**GODogs**](https://www.godogs.org.uk/) Project.

# Overview
This is an exploration of how imputable array-level data might be against the powerful reference dataset of the *Dog10k* dataset[^1].
In order to test this we take a set of variants deteted from the sequencing of 676 canines[^2] which we will refer to as the *Ostrander* dataset.
Each VCF file from the *Ostrander* dataset is downsampled such that new VCFs containing a restricted number of variants present on the commercial [wisdom array](https://www.wisdompanel.com/en-gb)[^3].
The idea here is to downsample *Ostrander* to the level of the wisdom array and then impute back up using *Dog10k* as the reference panel.

We can then explore how good the imputation process was for data obtained from a relatively inexpensive commercial array platform.

## HPC Setup
These analyses were run on the *CGS* production cluster (`sinfo -l -N -r `): 

```
NODELIST   NODES PARTITION       STATE CPUS    S:C:T MEMORY TMP_DISK WEIGHT AVAIL_FE REASON              
planex0        1  PlanEx2*       mixed 48     2:12:2 500000        0      1   (null) none                
planex1        1  PlanEx2*       mixed 48     2:12:2 500000        0      1   (null) none                
planex2        1  PlanEx2*       mixed 48     2:12:2 500000        0      1   (null) none
planex3        1  PlanEx2*       mixed 48     2:12:2 500000        0      1   (null) none                
planex4        1  PlanEx2*       mixed 48     2:12:2 500000        0      1   (null) none                
planex5        1  PlanEx2*       mixed 48     2:12:2 500000        0      1   (null) none
planex6        1  PlanEx2*       mixed 48     2:12:2 500000        0      1   (null) none
```

This cluster has 336 CPUs across 7 nodes, each with 500Gb RAM and is controlled using the *slurm* workload manager. For this analysis nodes 3-6 were out of action.
All analysis was hence carried out on nodes 0-2.

## Data Sources

### Ostrander Dataset
The Ostrander dataset[^2] was provided by the Raffan Lab as individual chromosome specific VCF files (One for each studied chromosome).
The Ostrander data from Broad is based on *canFam3*. We will need to liftover the SNV coordinates to *canFam4* to match our reference.

<details>
<summary>List of Input Files</summary>
  
### Ostrander Files
  
```
broad-chr1.vcf.gz
broad-chr2.vcf.gz
broad-chr3.vcf.gz
broad-chr4.vcf.gz
broad-chr5.vcf.gz
broad-chr6.vcf.gz
broad-chr7.vcf.gz
broad-chr8.vcf.gz
broad-chr9.vcf.gz
broad-chr10.vcf.gz
broad-chr11.vcf.gz
broad-chr12.vcf.gz
broad-chr13.vcf.gz
broad-chr14.vcf.gz
broad-chr15.vcf.gz
broad-chr16.vcf.gz
broad-chr17.vcf.gz
broad-chr18.vcf.gz
broad-chr19.vcf.gz
broad-chr20.vcf.gz
broad-chr21.vcf.gz
broad-chr22.vcf.gz
broad-chr23.vcf.gz
broad-chr24.vcf.gz
broad-chr25.vcf.gz
broad-chr26.vcf.gz
broad-chr27.vcf.gz
broad-chr28.vcf.gz
broad-chr29.vcf.gz
broad-chr30.vcf.gz
broad-chr31.vcf.gz
broad-chr32.vcf.gz
broad-chr33.vcf.gz
broad-chr34.vcf.gz
broad-chr35.vcf.gz
broad-chr36.vcf.gz
broad-chr37.vcf.gz
broad-chr38.vcf.gz
```

</details>

### Dog10K Dataset
The *Dog10k* data[^1] was downloaded as *plink* compatible files from their [website](https://dog10k.kiz.ac.cn/Home/Download). 
This dataset is based on a new reference genome *CanFam4*. 
The following files were obtained (1.4Tb of data): 

|Size | Filename |
|-----|----------|
|240K | dog10K-alignment-sample-table.2022-02-23-v7.txt|
|16K  | dog10k_SNPs.log|
|782M | dog10k_SNPs.map|
|48K | dog10k_SNPs.nosex|
|393G | dog10k_SNPs.ped|
|249G | dog10k_SNPs.ped.1|
|25G | dog10k.SNPs.plink.bed|
|984M | dog10k.SNPs.plink.bim|
|64K | dog10k.SNPs.plink.fam|
|48K | dog10k.SNPs.plink.nosex|

### Wisdom Array Postions
Individual coordinates from the commercial wisdom array were also provided [markers_v5.csv](wisdom/markers_v5.csv).
These should be *canFam3* coordinates.

### Canine Karyotype
For the imputation later we want the official length of each canine chromosome from *canFam4*. We obtain this from [Ensembl](https://ftp.ensembl.org/pub/release-114/tsv/canis_lupus_familiaris/).
The file used is [Canis_lupus_familiaris.ROS_Cfam_1.0.114.karyotype.tsv](etc/Canis_lupus_familiaris.ROS_Cfam_1.0.114.karyotype.tsv).

### LiftOver Data

We download the liftover chain from UCSC to convert *canFam3* to *canFam4* later to make our datasets compatible.
We also need the raw FASTA assemblies for both and we prepare the data as follows for later:

```
wget https://hgdownload.soe.ucsc.edu/goldenPath/canFam3/liftOver/canFam3ToCanFam4.over.chain.gz
wget https://hgdownload.soe.ucsc.edu/goldenPath/canFam3/bigZips/canFam3.fa.gz
wget https://hgdownload.soe.ucsc.edu/goldenPath/canFam4/bigZips/canFam4.fa.gz
gunzip canFam3.fa.gz canFam4.fa.gz
samtools faidx canFam3.fa
samtools faidx canFam4.fa
```

### Data Overview 

| Panel | SNVs | Individuals | Genome Build |
|-------|------|-------------|--------------|
|Dog10K[^1] | 51M  | 1987 | CanFam 4            |
|Ostrander[^2]| 91M | 676 | CanFam 3.1          |
|Wisdom[^3]|95,165 |  | CanFam 3.1              |


## Software Used

* [plink](https://www.cog-genomics.org/plink/) - PLINK v1.90b6.24 64-bit (6 Jun 2021)
* [plink2](https://www.cog-genomics.org/plink/2.0/) - PLINK v2.0.0-a.6.16LM AVX2 Intel (9 Jun 2025)
* [bcftools](https://github.com/samtools/bcftools) - 1.19-46-gc63329bb
* [samtools](https://github.com/samtools/samtools) - 1.19.2-7-g0c03cbd
* [score](https://github.com/freeseek/score) - v1 2024[^5]
* [shapeit](https://mathgen.stats.ox.ac.uk/genetics_software/shapeit/shapeit.html) - Version : v2.r904 
* [impute2](https://mathgen.stats.ox.ac.uk/impute/impute_v2.html) - version 2.3.2

---
# Imputation Pipeline
---

We first prepare the *Ostrander* dataset which is already split into Chromosomes. 

### 1. Ostrander Preparation
For each Chromosome XX (Where XX is 1-38) grab its VCF file we do the following:

* Downsample to the SNVs from that chromosome based on the wisdom positions (see above) using [downsample.pl](scripts/downsample.pl).
* Liftover the coordinates from *CanFam3* to *CanFam4* to match our reference.

#### 1a. Downsampling
* Downsampling is run by a custom perl script [downsample.pl](scripts/downsample.pl).
  * Invocation: `./downsample.pl broad-chr*.vcf.gz`
  * Takes a VCF file for a particular chromosome and removes positions not present in the [markers_v5.csv](wisdom/markers_v5.csv) file.
  * Produces a new file with `.downsampled.vcf` extension and gzips it to `.downsampled.vcf.gz`.


<details>
<summary>List of Downsampled Files Produced</summary>
  
### Ostrander Files
  
```
broad-chr1.downsampled.vcf.gz
broad-chr2.downsampled.vcf.gz
broad-chr3.downsampled.vcf.gz
broad-chr4.downsampled.vcf.gz
broad-chr5.downsampled.vcf.gz
broad-chr6.downsampled.vcf.gz
broad-chr7.downsampled.vcf.gz
broad-chr8.downsampled.vcf.gz
broad-chr9.downsampled.vcf.gz
broad-chr10.downsampled.vcf.gz
broad-chr11.downsampled.vcf.gz
broad-chr12.downsampled.vcf.gz
broad-chr13.downsampled.vcf.gz
broad-chr14.downsampled.vcf.gz
broad-chr15.downsampled.vcf.gz
broad-chr16.downsampled.vcf.gz
broad-chr17.downsampled.vcf.gz
broad-chr18.downsampled.vcf.gz
broad-chr19.downsampled.vcf.gz
broad-chr20.downsampled.vcf.gz
broad-chr21.downsampled.vcf.gz
broad-chr22.downsampled.vcf.gz
broad-chr23.downsampled.vcf.gz
broad-chr24.downsampled.vcf.gz
broad-chr25.downsampled.vcf.gz
broad-chr26.downsampled.vcf.gz
broad-chr27.downsampled.vcf.gz
broad-chr28.downsampled.vcf.gz
broad-chr29.downsampled.vcf.gz
broad-chr30.downsampled.vcf.gz
broad-chr31.downsampled.vcf.gz
broad-chr32.downsampled.vcf.gz
broad-chr33.downsampled.vcf.gz
broad-chr34.downsampled.vcf.gz
broad-chr35.downsampled.vcf.gz
broad-chr36.downsampled.vcf.gz
broad-chr37.downsampled.vcf.gz
broad-chr38.downsampled.vcf.gz
```

</details>

#### 1b. Liftover

We perform the Liftover for the downsampled VCF files from *canFam3* to *canFam4* using the UCSC liftover path (described above).
For this we will process each VCF, sanitise it (the chrXX ids in particular), then use the *bcftools* liftover plugin from the *score* toolkit[^5].
We provide the liftover chain path `canFam3ToCanFam4.over.chain.gz` and the two reference FASTA files (`canFam3.fa` and `canFam4.fa`).

The command works something like:
* `bcftools +liftover --no-version -Ou /tmp/broad-chrXX.downsampled.vcf -- -c ../liftover/canFam3ToCanFam4.over.chain.gz -s ../liftover/canFam3.fa -f ../liftover/canFam4.fa | bcftools sort -Oz -o /tmp/broad-chrXX.downsampled.cf4.vcf.gz`

This is run using a custom perl script [liftover.pl](scripts/liftover.pl) that operates as follows:
* Script invovation: `liftover.pl *.downsampled.vcf.gz`
* Runtime: *Approximately 15 secs per chromosome, 1cpu*
* For each downsampled VCF:
  * Decompress and fix the CHR naming nomenclature, create a temporary VCF in /tmp
  * Run the *bcftools* with liftover plugin against the two references using the UCSC liftover path from *canFam3* to *canFam4* to a new cf4.vcf.gz file.
  * Decompress and then fix back the CHR naming nomenclature to a new VCF.gz file that's copied back to the path from /tmp
    

<details>
<summary>List of Downsampled LiftOver Files Produced</summary>
  
### Ostrander Liftover Files canFam4
  
```
broad-chr1.downsampled.cf4.vcf.gz
broad-chr2.downsampled.cf4.vcf.gz
broad-chr3.downsampled.cf4.vcf.gz
broad-chr4.downsampled.cf4.vcf.gz
broad-chr5.downsampled.cf4.vcf.gz
broad-chr6.downsampled.cf4.vcf.gz
broad-chr7.downsampled.cf4.vcf.gz
broad-chr8.downsampled.cf4.vcf.gz
broad-chr9.downsampled.cf4.vcf.gz
broad-chr10.downsampled.cf4.vcf.gz
broad-chr11.downsampled.cf4.vcf.gz
broad-chr12.downsampled.cf4.vcf.gz
broad-chr13.downsampled.cf4.vcf.gz
broad-chr14.downsampled.cf4.vcf.gz
broad-chr15.downsampled.cf4.vcf.gz
broad-chr16.downsampled.cf4.vcf.gz
broad-chr17.downsampled.cf4.vcf.gz
broad-chr18.downsampled.cf4.vcf.gz
broad-chr19.downsampled.cf4.vcf.gz
broad-chr20.downsampled.cf4.vcf.gz
broad-chr21.downsampled.cf4.vcf.gz
broad-chr22.downsampled.cf4.vcf.gz
broad-chr23.downsampled.cf4.vcf.gz
broad-chr24.downsampled.cf4.vcf.gz
broad-chr25.downsampled.cf4.vcf.gz
broad-chr26.downsampled.cf4.vcf.gz
broad-chr27.downsampled.cf4.vcf.gz
broad-chr28.downsampled.cf4.vcf.gz
broad-chr29.downsampled.cf4.vcf.gz
broad-chr30.downsampled.cf4.vcf.gz
broad-chr31.downsampled.cf4.vcf.gz
broad-chr32.downsampled.cf4.vcf.gz
broad-chr33.downsampled.cf4.vcf.gz
broad-chr34.downsampled.cf4.vcf.gz
broad-chr35.downsampled.cf4.vcf.gz
broad-chr36.downsampled.cf4.vcf.gz
broad-chr37.downsampled.cf4.vcf.gz
broad-chr38.downsampled.cf4.vcf.gz
```

</details>

#### 1c. PLINK conversion and filtering
* Prepare a PLINK dataset for each downsampled Chromosome (XX where XX is Chr1 to Chr38).
  * `plink --const-fid 0 --vcf ../broad-chrXX.downsampled.cf4.vcf.gz --out broad_plink_chrXX --dog`
* Fix the Name column as it still contains CHR:POS names from *cf3.1*:
  * `perl -lane 'next if($F[0] eq "0"); print $F[1]."\t".$F[0].":".$F[3];' broad_plink_chrXX.bim > broad_plink_chrXX.names`
  * `plink --update-name broad_plink_chrXX.names --make-bed --bfile broad_plink_chrXX --out broad_plink_chrXX.1 --dog`
* Find and store ambiguous SNPs from the new plink file
  * `perl -lane 'if(($F[4] eq "A" && $F[5] eq "T") || ($F[4] eq "T" && $F[5] eq "A") || ($F[4] eq "C" && $F[5] eq "G") || ($F[4] eq "G" && $F[5] eq "C")){ print $F[1] }' broad_plink_chrXX.1.bim > ambiguous.snps.chrXX`
* Exclude the ambiguous SNPs in a new plink file
  * `plink --bfile broad_plink_chrXX.1 --exclude ambiguous.snps.chrXX --make-bed --out broad_plink_chrXX.2 --dog`
* Continue with filtering - We adjusted the Hard-Weinberg P-value filter as too many SNVs were being excluded given the size of the dataset.
  * `plink --bfile broad_plink_chrXX.2 --geno 0.03 --make-bed --out broad_plink_chrXX.geno --dog`
  * `plink --bfile broad_plink_chrXX.geno --mind 0.1 --make-bed --out broad_plink_chrXX.mind --dog`
  * `plink --bfile broad_plink_chrXX.mind --maf 0.01 --make-bed --out broad_plink_chrXX.maf --dog`
  * `plink --bfile broad_plink_chrXX.maf --hwe 0.000000000000000000005 --make-bed --out broad_plink_chrXX_gwas --dog`


### HPC Details:
* Two scripts are used for this:
  * [run_chr.sh](scripts/run_chr.sh) - A generic script that runs a downsampled chromosome through the plink commands above.
  * [run_chr.pl](scripts/run_chr.pl) - A perl script that generates an *sbatch* script for each chromosome that will launch *run_chr.sh* via *sbatch* on our HPC.
* **Runtime - Downsampling, 5-10 mins per chromosome, PLINK: approximately 30-60secs per chromosome**.
  * sbatch parameters used:
    * `--partition=PlanEx2` (use the CGS cluster).
    * `--mem=2G` (2Gb RAM needed)
    * `--ntasks=1` (1 cpu per job)
    * `-o XXXX,out` (output log)
    * `-e YYYY.err` (error log)
 

<details>
<summary>This should produce the following files for each downsampled chromosome (example Chr10 below)</summary>

```
broad_plink_chr10.fam
broad_plink_chr10.bim
broad_plink_chr10.bed
broad_plink_chr10.log
broad_plink_chr10.2.nosex
broad_plink_chr10.2.bed
broad_plink_chr10.2.fam
broad_plink_chr10.2.bim
broad_plink_chr10.2.log
broad_plink_chr10.geno.nosex
broad_plink_chr10.geno.bed
broad_plink_chr10.geno.fam
broad_plink_chr10.geno.bim
broad_plink_chr10.geno.log
broad_plink_chr10.mind.nosex
broad_plink_chr10.mind.irem
broad_plink_chr10.mind.bed
broad_plink_chr10.mind.fam
broad_plink_chr10.mind.bim
broad_plink_chr10.mind.log
broad_plink_chr10.maf.nosex
broad_plink_chr10.maf.bed
broad_plink_chr10.maf.fam
broad_plink_chr10.maf.bim
broad_plink_chr10.maf.log
broad_plink_chr10_gwas.nosex
broad_plink_chr10_gwas.bed
broad_plink_chr10_gwas.fam
broad_plink_chr10_gwas.bim
broad_plink_chr10_gwas.log
```

</details>

<details>
<summary>Example Output from Conversion and Filtering Step (chr10)</summary>

```
../downsampled_liftover/broad-chr10.downsampled.cf4.vcf.gz
10
PLINK v1.90b6.24 64-bit (6 Jun 2021)           www.cog-genomics.org/plink/1.9/
(C) 2005-2021 Shaun Purcell, Christopher Chang   GNU General Public License v3
Logging to broad_plink_chr10.log.
Options in effect:
  --allow-extra-chr
  --chr 10
  --const-fid 0
  --dog
  --out broad_plink_chr10
  --vcf ../downsampled_liftover/broad-chr10.downsampled.cf4.vcf.gz

515857 MB RAM detected; reserving 257928 MB for main workspace.
--vcf: broad_plink_chr10.bed + broad_plink_chr10.bim + broad_plink_chr10.fam
written.
PLINK v1.90b6.24 64-bit (6 Jun 2021)           www.cog-genomics.org/plink/1.9/
(C) 2005-2021 Shaun Purcell, Christopher Chang   GNU General Public License v3
Logging to broad_plink_chr10.1.log.
Options in effect:
  --bfile broad_plink_chr10
  --dog
  --make-bed
  --out broad_plink_chr10.1
  --update-name broad_plink_chr10.names

515857 MB RAM detected; reserving 257928 MB for main workspace.
2475 variants loaded from .bim file.
676 dogs (0 males, 0 females, 676 ambiguous) loaded from .fam.
Ambiguous sex IDs written to broad_plink_chr10.1.nosex .
--update-name: 2475 values updated.
Using 1 thread (no multithreaded calculations invoked).
Before main variant filters, 676 founders and 0 nonfounders present.
Calculating allele frequencies... done.
Total genotyping rate is 0.97575.
2475 variants and 676 dogs pass filters and QC.
Note: No phenotypes present.
--make-bed to broad_plink_chr10.1.bed + broad_plink_chr10.1.bim +
broad_plink_chr10.1.fam ... done.
PLINK v1.90b6.24 64-bit (6 Jun 2021)           www.cog-genomics.org/plink/1.9/
(C) 2005-2021 Shaun Purcell, Christopher Chang   GNU General Public License v3
Logging to broad_plink_chr10.2.log.
Options in effect:
  --bfile broad_plink_chr10.1
  --dog
  --exclude broad_plink_ambiguous.snps.chr10
  --make-bed
  --out broad_plink_chr10.2

515857 MB RAM detected; reserving 257928 MB for main workspace.
2475 variants loaded from .bim file.
676 dogs (0 males, 0 females, 676 ambiguous) loaded from .fam.
Ambiguous sex IDs written to broad_plink_chr10.2.nosex .
--exclude: 2314 variants remaining.
Using 1 thread (no multithreaded calculations invoked).
Before main variant filters, 676 founders and 0 nonfounders present.
Calculating allele frequencies... done.
Total genotyping rate is 0.975601.
2314 variants and 676 dogs pass filters and QC.
Note: No phenotypes present.
--make-bed to broad_plink_chr10.2.bed + broad_plink_chr10.2.bim +
broad_plink_chr10.2.fam ... done.
PLINK v1.90b6.24 64-bit (6 Jun 2021)           www.cog-genomics.org/plink/1.9/
(C) 2005-2021 Shaun Purcell, Christopher Chang   GNU General Public License v3
Logging to broad_plink_chr10.geno.log.
Options in effect:
  --bfile broad_plink_chr10.2
  --dog
  --geno 0.03
  --make-bed
  --out broad_plink_chr10.geno

515857 MB RAM detected; reserving 257928 MB for main workspace.
2314 variants loaded from .bim file.
676 dogs (0 males, 0 females, 676 ambiguous) loaded from .fam.
Ambiguous sex IDs written to broad_plink_chr10.geno.nosex .
Using 1 thread (no multithreaded calculations invoked).
Before main variant filters, 676 founders and 0 nonfounders present.
Calculating allele frequencies... done.
Total genotyping rate is 0.975601.
444 variants removed due to missing genotype data (--geno).
1870 variants and 676 dogs pass filters and QC.
Note: No phenotypes present.
--make-bed to broad_plink_chr10.geno.bed + broad_plink_chr10.geno.bim +
broad_plink_chr10.geno.fam ... done.
PLINK v1.90b6.24 64-bit (6 Jun 2021)           www.cog-genomics.org/plink/1.9/
(C) 2005-2021 Shaun Purcell, Christopher Chang   GNU General Public License v3
Logging to broad_plink_chr10.mind.log.
Options in effect:
  --bfile broad_plink_chr10.geno
  --dog
  --make-bed
  --mind 0.1
  --out broad_plink_chr10.mind

515857 MB RAM detected; reserving 257928 MB for main workspace.
1870 variants loaded from .bim file.
676 dogs (0 males, 0 females, 676 ambiguous) loaded from .fam.
Ambiguous sex IDs written to broad_plink_chr10.mind.nosex .
6 dogs removed due to missing genotype data (--mind).
IDs written to broad_plink_chr10.mind.irem .
Using 1 thread (no multithreaded calculations invoked).
Before main variant filters, 670 founders and 0 nonfounders present.
Calculating allele frequencies... done.
Total genotyping rate in remaining samples is 0.988929.
1870 variants and 670 dogs pass filters and QC.
Note: No phenotypes present.
--make-bed to broad_plink_chr10.mind.bed + broad_plink_chr10.mind.bim +
broad_plink_chr10.mind.fam ... done.
PLINK v1.90b6.24 64-bit (6 Jun 2021)           www.cog-genomics.org/plink/1.9/
(C) 2005-2021 Shaun Purcell, Christopher Chang   GNU General Public License v3
Logging to broad_plink_chr10.maf.log.
Options in effect:
  --bfile broad_plink_chr10.mind
  --dog
  --maf 0.01
  --make-bed
  --out broad_plink_chr10.maf

515857 MB RAM detected; reserving 257928 MB for main workspace.
1870 variants loaded from .bim file.
670 dogs (0 males, 0 females, 670 ambiguous) loaded from .fam.
Ambiguous sex IDs written to broad_plink_chr10.maf.nosex .
Using 1 thread (no multithreaded calculations invoked).
Before main variant filters, 670 founders and 0 nonfounders present.
Calculating allele frequencies... done.
Total genotyping rate is 0.988929.
11 variants removed due to minor allele threshold(s)
(--maf/--max-maf/--mac/--max-mac).
1859 variants and 670 dogs pass filters and QC.
Note: No phenotypes present.
--make-bed to broad_plink_chr10.maf.bed + broad_plink_chr10.maf.bim +
broad_plink_chr10.maf.fam ... done.
PLINK v1.90b6.24 64-bit (6 Jun 2021)           www.cog-genomics.org/plink/1.9/
(C) 2005-2021 Shaun Purcell, Christopher Chang   GNU General Public License v3
Logging to broad_plink_chr10_gwas.log.
Options in effect:
  --bfile broad_plink_chr10.maf
  --dog
  --hwe 0.000000000000000000005
  --make-bed
  --out broad_plink_chr10_gwas

515857 MB RAM detected; reserving 257928 MB for main workspace.
1859 variants loaded from .bim file.
670 dogs (0 males, 0 females, 670 ambiguous) loaded from .fam.
Ambiguous sex IDs written to broad_plink_chr10_gwas.nosex .
Using 1 thread (no multithreaded calculations invoked).
Before main variant filters, 670 founders and 0 nonfounders present.
Calculating allele frequencies... done.
Total genotyping rate is 0.988925.
--hwe: 153 variants removed due to Hardy-Weinberg exact test.
1706 variants and 670 dogs pass filters and QC.
Note: No phenotypes present.
--make-bed to broad_plink_chr10_gwas.bed + broad_plink_chr10_gwas.bim +
broad_plink_chr10_gwas.fam ... done.
```
</details>

---

### Filtering Statistics

From our original wisdom panel downsampled VCF files we can see how many SNVs we lose through liftover and then *PLINK* filtering.
This could be more or less stringent, but this seems sensible for our goals here.

#### SNV Table
|Chromosome|Downsampled VCF| Liftover VCF | Final Filtered|
|----------|-----------|--------------|---------------|
|chr1|4326|4324|3134|
|chr2|2946|2945|2133|
|chr3|3409|3408|2514|
|chr4|3232|3232|2489|
|chr5|3266|3264|2316|
|chr6|2845|2842|2110|
|chr7|2984|2984|2338|
|chr8|2586|2586|2031|
|chr9|2138|2136|1363|
|chr10|2475|2475|1706|
|chr11|2443|2441|1855|
|chr12|2834|2832|2161|
|chr13|2377|2375|1847|
|chr14|2149|2148|1760|
|chr15|2225|2225|1608|
|chr16|2129|2126|1643|
|chr17|2352|2350|1815|
|chr18|2045|2044|1395|
|chr19|1962|1961|1594|
|chr20|2084|2083|1491|
|chr21|1868|1865|1499|
|chr22|2254|2254|1822|
|chr23|1991|1991|1635|
|chr24|1847|1846|1301|
|chr25|1964|1963|1489|
|chr26|1472|1472|1017|
|chr27|1744|1744|1407|
|chr28|1567|1567|1109|
|chr29|1514|1514|1259|
|chr30|1563|1563|1191|
|chr31|1426|1424|1085|
|chr32|1429|1428|1185|
|chr33|1181|1181|958|
|chr34|1697|1695|1349|
|chr35|1105|1104|858|
|chr36|1118|1117|925|
|chr37|1140|1137|880|
|chr38|850|847|645|
|**SUM**|**80537**|**80493**|**60917**|


---

### 2. Ostrander Phasing

For each Chromosome XX we now perform the phasing analysis using *shapeit*, again we use a HPC script (below) and we also need a genetic map.
The [genetic map used](https://github.com/cflerin/dog_recombination)[^4] has a map of each Chromosome, this will be used for all phasing and imputation.
  * The maps used are in the [maps folder](./maps).
  * `shapeit -M ../../maps/chrXX.cf3.1_map.txt -B broad_plink_chrXX_gwas -O broad_plink_chrXX_gwas.phased  -T 4 --window 2 --effective-size 200 --force`
  * `shapeit -convert --input-haps broad_plink_chrXX_gwas.phased --output-ref broad_plink_chrXX_gwas.phased.impute`

### HPC Details:
* **Runtime - Approximately 10 mins per chromosome.**
* One script is used for this:
  * [shape_it_gwas.pl](scripts/shape_it_gwas.pl) - A perl script that generates an *sbatch* script for each chromosome that will launch via *sbatch* on our HPC.
  * run without parameters it will print the commands and scripts without submitting.
  * To actually submit the jobs run with the `--runnit` flag, e.g. `./shape_it_gwas.pl --runnit`
  * sbatch parameters used:
    * `--partition=PlanEx2` (use the CGS cluster).
    * `--mem=15G` (15Gb RAM needed)
    * `--ntasks=1` (1 cpu per job)
    * `--cpus-per-task=4` (4 threads per job)
    * `-o XXXX,out` (output log)
    * `-e YYYY.err` (error log)

<details>
<summary>This should produce the following new files for each downsampled chromosome (example Chr10 below)</summary>

```
broad_plink_chr10_gwas.phased.sample
broad_plink_chr10_gwas.phased.haps
broad_plink_chr10_gwas.phased.impute.samples
broad_plink_chr10_gwas.phased.impute.legend
broad_plink_chr10_gwas.phased.impute.haplotypes
```

</details>

<details>
<summary>Example Phasing Output (chr10)</summary>

```

Segmented HAPlotype Estimation & Imputation Tool
  * Authors : Olivier Delaneau, Jared O'Connell, Jean-François Zagury, Jonathan Marchini
  * Contact : send an email to the OXSTATGEN mail list https://www.jiscmail.ac.uk/cgi-bin/webadmin?A0=OXSTATGEN
  * Webpage : https://mathgen.stats.ox.ac.uk/shapeit
  * Version : v2.r904
  * Date    : 13/06/2025 21:15:23
  * LOGfile : [shapeit_13062025_21h15m23s_31cdafc1-d539-4187-9ad0-505c15dbf9d2.log]

MODE -phase : PHASING GENOTYPE DATA
  * Autosome (chr1 ... chr22)
  * Window-based model (SHAPEIT v2)
  * MCMC iteration

Parameters :
  * Seed : 1749845723
  * Parallelisation: 4 threads
  * Ref allele is NOT aligned on the reference genome
  * MCMC: 35 iterations [7 B + 1 runs of 8 P + 20 M]
  * Model: 100 states per window [100 H + 0 PM + 0 R + 0 COV ] / Windows of ~2.0 Mb / Ne = 200

Reading site list in [broad_plink_chr10_gwas.bim]
  * 1706 sites included

Reading sample list in [broad_plink_chr10_gwas.fam]
  * 670 samples included
  * 670 unrelateds / 0 duos / 0 trios in 670 different families

Reading genotypes in [broad_plink_chr10_gwas.bed]
  * Plink binary file SNP-major mode

Reading genetic map in [../../maps/chr10.cf3.1_map.txt]
  * 345 genetic positions found
  * #set=0 / #interpolated=1706
  * Physical map [0.45 Mb -> 70.54 Mb] / Genetic map [0.00 cM -> 63.83 cM]

Checking missingness and MAF...
  * 33 individuals with high rates of missing data (>5%)

Building graphs [670/670]
  * 670 graphs / 114446 segments / ~9 SNPs per segment / 4407501 transitions
  * 0 haploids / 670 unrelateds / 0 duos / 0 trios
  * 1340 founder haplotypes

Sampling haplotypes [670/670]

Burn-in iteration [1/7] [670/670]

Burn-in iteration [2/7] [670/670]

Burn-in iteration [3/7] [670/670]

Burn-in iteration [4/7] [670/670]

Burn-in iteration [5/7] [670/670]

Burn-in iteration [6/7] [670/670]

Burn-in iteration [7/7] [670/670]

Pruning iteration [1/8] [670/670]

Pruning iteration [2/8] [670/670]

Pruning iteration [3/8] [670/670]

Pruning iteration [4/8] [670/670]

Pruning iteration [5/8] [670/670]

Pruning iteration [6/8] [670/670]

Pruning iteration [7/8] [670/670]

Pruning iteration [8/8] [670/670]

Pruning graphs [670/670]
  * 670 graphs / 57790 segments / ~19 SNPs per segment / 1368412 transitions
  * 0 haploids / 670 unrelateds / 0 duos / 0 trios
  * 1340 founder haplotypes

Main iteration [1/20] [670/670]

Main iteration [2/20] [670/670]

Main iteration [3/20] [670/670]

Main iteration [4/20] [670/670]

Main iteration [5/20] [670/670]

Main iteration [6/20] [670/670]

Main iteration [7/20] [670/670]

Main iteration [8/20] [670/670]

Main iteration [9/20] [670/670]

Main iteration [10/20] [670/670]

Main iteration [11/20] [670/670]

Main iteration [12/20] [670/670]

Main iteration [13/20] [670/670]

Main iteration [14/20] [670/670]

Main iteration [15/20] [670/670]

Main iteration [16/20] [670/670]

Main iteration [17/20] [670/670]

Main iteration [18/20] [670/670]

Main iteration [19/20] [670/670]

Main iteration [20/20] [670/670]

Normalising graphs [670/670]

Solving haplotypes [670/670]

Writing sample list in [broad_plink_chr10_gwas.phased.sample]

Writing site list and haplotypes in [broad_plink_chr10_gwas.phased.haps]

Running time: 687 seconds

Segmented HAPlotype Estimation & Imputation Tool
  * Authors : Olivier Delaneau, Jared O'Connell, Jean-François Zagury, Jonathan Marchini
  * Contact : send an email to the OXSTATGEN mail list https://www.jiscmail.ac.uk/cgi-bin/webadmin?A0=OXSTATGEN
  * Webpage : https://mathgen.stats.ox.ac.uk/shapeit
  * Version : v2.r904
  * Date    : 13/06/2025 21:26:50
  * LOGfile : [shapeit_13062025_21h26m50s_efa2ce3f-a959-4fac-909f-0e9ecb760fde.log]

MODE -convert : CONVERTING/SUBSETTING SHAPEIT HAPLOTYPES

Parameters :
  * Seed : 1749846410
  * Parallelisation: 1 threads
  * Ref allele is NOT aligned on the reference genome

Reading sample list in [broad_plink_chr10_gwas.phased.sample]
  * 670 individuals included

Reading site list and haplotypes in [broad_plink_chr10_gwas.phased.haps]
  * 1706 sites included

Mapping sites and samples of the haplotype set
  * 1706 sites found
  * 0 sites not found
  * 1340 founder haplotypes

Writing sample list in [broad_plink_chr10_gwas.phased.impute.samples]

Writing site list in [broad_plink_chr10_gwas.phased.impute.legend]

Writing haplotypes in [broad_plink_chr10_gwas.phased.impute.haplotypes]

Running time: 1 seconds

```
</details>


---

### 3. Dog10K Preparation

We will take the giant PLINK formatted datafile `dog10k.SNPs.plink` from the Dog10K dataset and split it first into individual chromosome (XX is the Chr number), PLINK files.
We will also do some simple filtering and create a new *refpanel* plink object for each.
  * `plink --bfile ../dog10k.SNPs.plink --chr XX --make-bed --out dog10k_plink_chrXX --dog`
  * `plink --bfile dog10k_plink_chrXX --maf 0.01 --mind 0.1 --geno 0.03 --make-bed --out dog10k_plink_chrXX.refpanel --dog`


Finally we want to fix the names and generate Minor Allele Frequencies for later.
  * `plink2 --bfile dog10k_plink_chrXX.refpanel --set-all-var-ids @:# --make-bed --out dog10k_plink_chrXX.refpanel.names`
  * `plink2 --bfile dog10k_plink_chrXX.refpanel.names --freq --out dog10k_plink_chrXX.refpanel.names.maf`

This generates files like this for Chr31: `dog10k_plink_chr31.refpanel.names.maf.afreq`.
```
#CHROM  ID      REF     ALT     PROVISIONAL_REF?        ALT_FREQS       OBS_CT
31      31:740  C       A       Y       0.0543807       3972
31      31:822  T       A       Y       0.166667        3972
31      31:925  G       A       Y       0.0468278       3972
31      31:945  A       G       Y       0.167674        3972
31      31:1024 A       G       Y       0.0370091       3972
31      31:1186 C       T       Y       0.0239174       3972
31      31:1192 G       T       Y       0.0239174       3972
31      31:1215 G       T       Y       0.0475831       3972
31      31:1343 T       C       Y       0.0460957       3970
31      31:1404 T       C       Y       0.0306995       3974
31      31:1529 G       A       Y       0.0362355       3974
31      31:1630 C       T       Y       0.0405133       3974
31      31:1647 G       T       Y       0.0410166       3974
```
We'll use this later for the concordance analysis.


### HPC Details:
* **Runtime - Approximately 10-20 mins per chromosome.**
* One script is used for this:
  * [process_chr.pl](scripts/process_chr.pl) - A perl script that generates an *sbatch* script for each chromosome that will launch via *sbatch* on our HPC.
  * run without parameters it will print the commands and scripts without submitting.
  * To actually submit the jobs run with the `--runnit` flag, e.g. `./process_chr.pl --runnit`
  * sbatch parameters used:
    * `--partition=PlanEx2` (use the CGS cluster).
    * `--mem=10G` (10Gb RAM needed)
    * `--ntasks=1` (1 cpu per job)
    * `--cpus-per-task=1` (1 thread per job)
    * `-o XXXX,out` (output log)
    * `-e YYYY.err` (error log)

<details>
<summary>This should produce the following new files for each Dog10K chromosome (example Chr10 below)</summary>

```
dog10k_plink_chr10.nosex
dog10k_plink_chr10.bed
dog10k_plink_chr10.fam
dog10k_plink_chr10.bim
dog10k_plink_chr10.log
dog10k_plink_chr10.refpanel.nosex
dog10k_plink_chr10.refpanel.bed
dog10k_plink_chr10.refpanel.fam
dog10k_plink_chr10.refpanel.bim
dog10k_plink_chr10.refpanel.log
```

</details>

---

### 4. Dog10K Phasing

For each Dog10K Chromosome XX we now perform the phasing analysis using *shapeit*, again we use a HPC script (below) and we also need a genetic map.
The [genetic map used](https://github.com/cflerin/dog_recombination)[^4] has a map of each Chromosome, this will be used for all phasing and imputation.
  * The maps used are in the [maps folder](./maps).
  * `shapeit -M ../../maps/chrXX.cf3.1_map.txt -B dog10k_plink_chrXX.refpanel -O dog10k_plink_chrXX.phased -T 20 --window 2 --effective-size 200 --force`
  * `shapeit -convert --input-haps dog10k_plink_chrXX.phased --output-ref dog10k_plink_chrXX.phased.impute`

### HPC Details:
* **Runtime - Approximately 14-36 hours per chromosome on 20 cpus each.**
* One script is used for this:
  * [shape_it.pl](scripts/shape_it.pl) - A perl script that generates an *sbatch* script for each chromosome that will launch via *sbatch* on our HPC.
  * run without parameters it will print the commands and scripts without submitting.
  * To actually submit the jobs run with the `--runnit` flag, e.g. `./shape_it_gwas.pl --runnit`
  * sbatch parameters used:
    * `--partition=PlanEx2` (use the CGS cluster).
    * `--mem=60G` (60Gb RAM needed)
    * `--ntasks=1` (1 cpu per job)
    * `--cpus-per-task=20` (20 threads per job)
    * `-o XXXX,out` (output log)
    * `-e YYYY.err` (error log)

* A monitoring script [progress.pl](scripts/progress.pl) can update you on the status of these long running jobs:

<details>
  <summary>Example Monitoring Output of a run in progress is below:</summary>


|Chromosome	|File	|Size (MB)	|Status	|
|-----------|-----|-----------|-------|
|1	|dog10k_plink_chr1.refpanel.bed	|333M	|Finished: 35.38 (hours)	|
|2	|dog10k_plink_chr2.refpanel.bed	|230M	|Finished: 19.53 (hours)	|
|3	|dog10k_plink_chr3.refpanel.bed	|268M	|Finished: 25.98 (hours)	|
|4	|dog10k_plink_chr4.refpanel.bed	|248M	|Finished: 22.92 (hours)	|
|5	|dog10k_plink_chr5.refpanel.bed	|250M	|Finished: 23.02 (hours)	|
|6	|dog10k_plink_chr6.refpanel.bed	|211M	|Finished: 17.20 (hours)	|
|7	|dog10k_plink_chr7.refpanel.bed	|215M	|Finished: 17.94 (hours)	|
|8	|dog10k_plink_chr8.refpanel.bed	|214M	|Finished: 19.28 (hours)	|
|9	|dog10k_plink_chr9.refpanel.bed	|158M	|Finished: 11.86 (hours)	|
|10	|dog10k_plink_chr10.refpanel.bed	|190M	|Finished: 14.36 (hours)	|
|11	|dog10k_plink_chr11.refpanel.bed	|199M	|Finished: 15.95 (hours)	|
|12	|dog10k_plink_chr12.refpanel.bed	|216M	|Pruning graphs [1110/1987]	|
|13	|dog10k_plink_chr13.refpanel.bed	|190M	|Pruning graphs [1601/1987]	|
|14	|dog10k_plink_chr14.refpanel.bed	|165M	|Pruning graphs [1541/1987]	|
|15	|dog10k_plink_chr15.refpanel.bed	|173M	|Pruning graphs [551/1987]	|
|16	|dog10k_plink_chr16.refpanel.bed	|182M	|Pruning iteration [6/8] [280/1987]	|
|17	|dog10k_plink_chr17.refpanel.bed	|188M	|Pruning iteration [2/8] [1460/1987]	|
|18	|dog10k_plink_chr18.refpanel.bed	|173M	|Not Started	|
|19	|dog10k_plink_chr19.refpanel.bed	|167M	|Not Started	|
|20	|dog10k_plink_chr20.refpanel.bed	|154M	|Not Started	|
|21	|dog10k_plink_chr21.refpanel.bed	|164M	|Not Started	|
|22	|dog10k_plink_chr22.refpanel.bed	|172M	|Not Started	|
|23	|dog10k_plink_chr23.refpanel.bed	|155M	|Not Started	|
|24	|dog10k_plink_chr24.refpanel.bed	|145M	|Not Started	|
|25	|dog10k_plink_chr25.refpanel.bed	|159M	|Not Started	|
|26	|dog10k_plink_chr26.refpanel.bed	|145M	|Not Started	|
|27	|dog10k_plink_chr27.refpanel.bed	|145M	|Not Started	|
|28	|dog10k_plink_chr28.refpanel.bed	|132M	|Not Started	|
|29	|dog10k_plink_chr29.refpanel.bed	|135M	|Not Started	|
|30	|dog10k_plink_chr30.refpanel.bed	|118M	|Not Started	|
|31	|dog10k_plink_chr31.refpanel.bed	|146M	|Not Started	|
|32	|dog10k_plink_chr32.refpanel.bed	|137M	|Not Started	|
|33	|dog10k_plink_chr33.refpanel.bed	|99M	|Not Started	|
|34	|dog10k_plink_chr34.refpanel.bed	|143M	|Not Started	|
|35	|dog10k_plink_chr35.refpanel.bed	|112M	|Not Started	|
|36	|dog10k_plink_chr36.refpanel.bed	|91M	|Not Started	|
|37	|dog10k_plink_chr37.refpanel.bed	|95M	|Not Started	|
|38	|dog10k_plink_chr38.refpanel.bed	|98M	|Not Started	|


</details>


<details>
<summary>This should produce the following new files for each Dog10K chromosome (example Chr10 below)</summary>

```
dog10k_plink_chr10.phased.sample
dog10k_plink_chr10.phased.haps
dog10k_plink_chr10.phased.impute.samples
dog10k_plink_chr10.phased.impute.legend
dog10k_plink_chr10.impute.haplotypes
```

</details>

<details>
  <summary>Output Log Example (Chromosome 2)</summary>

```

Segmented HAPlotype Estimation & Imputation Tool
  * Authors : Olivier Delaneau, Jared O'Connell, Jean-François Zagury, Jonathan Marchini
  * Contact : send an email to the OXSTATGEN mail list https://www.jiscmail.ac.uk/cgi-bin/webadmin?A0=OXSTATGEN
  * Webpage : https://mathgen.stats.ox.ac.uk/shapeit
  * Version : v2.r904
  * Date    : 12/06/2025 12:11:08
  * LOGfile : [shapeit_12062025_12h11m08s_1afea299-f4be-4a66-b1a0-27bc00dbf0d2.log]

MODE -phase : PHASING GENOTYPE DATA
  * Autosome (chr1 ... chr22)
  * Window-based model (SHAPEIT v2)
  * MCMC iteration

Parameters :
  * Seed : 1749726668
  * Parallelisation: 20 threads
  * Ref allele is NOT aligned on the reference genome
  * MCMC: 35 iterations [7 B + 1 runs of 8 P + 20 M]
  * Model: 100 states per window [100 H + 0 PM + 0 R + 0 COV ] / Windows of ~2.0 Mb / Ne = 200

Reading site list in [dog10k_plink_chr2.refpanel.bim]
  * 483737 sites included

Reading sample list in [dog10k_plink_chr2.refpanel.fam]
  * 1987 samples included
  * 1987 unrelateds / 0 duos / 0 trios in 1987 different families

Reading genotypes in [dog10k_plink_chr2.refpanel.bed]
  * Plink binary file SNP-major mode

Reading genetic map in [../../maps/chr2.cf3.1_map.txt]
  * 419 genetic positions found
  * #set=3 / #interpolated=483734
  * Physical map [0.00 Mb -> 84.98 Mb] / Genetic map [0.00 cM -> 74.74 cM]

Checking missingness and MAF...

Building graphs [1987/1987]
  * 1987 graphs / 52990142 segments / ~18 SNPs per segment / 1764741174 transitions
  * 0 haploids / 1987 unrelateds / 0 duos / 0 trios
  * 3974 founder haplotypes

Sampling haplotypes [1987/1987]

Burn-in iteration [1/7] [1987/1987]

Burn-in iteration [2/7] [1987/1987]

Burn-in iteration [3/7] [1987/1987]

Burn-in iteration [4/7] [1987/1987]

Burn-in iteration [5/7] [1987/1987]

Burn-in iteration [6/7] [1987/1987]

Burn-in iteration [7/7] [1987/1987]

Pruning iteration [1/8] [1987/1987]

Pruning iteration [2/8] [1987/1987]

Pruning iteration [3/8] [1987/1987]

Pruning iteration [4/8] [1987/1987]

Pruning iteration [5/8] [1987/1987]

Pruning iteration [6/8] [1987/1987]

Pruning iteration [7/8] [1987/1987]

Pruning iteration [8/8] [1987/1987]

Pruning graphs [1987/1987]
  * 1987 graphs / 15733259 segments / ~61 SNPs per segment / 236024952 transitions
  * 0 haploids / 1987 unrelateds / 0 duos / 0 trios
  * 3974 founder haplotypes

Main iteration [1/20] [1987/1987]

Main iteration [2/20] [1987/1987]

Main iteration [3/20] [1987/1987]

Main iteration [4/20] [1987/1987]

Main iteration [5/20] [1987/1987]

Main iteration [6/20] [1987/1987]

Main iteration [7/20] [1987/1987]

Main iteration [8/20] [1987/1987]

Main iteration [9/20] [1987/1987]

Main iteration [10/20] [1987/1987]

Main iteration [11/20] [1987/1987]

Main iteration [12/20] [1987/1987]

Main iteration [13/20] [1987/1987]

Main iteration [14/20] [1987/1987]

Main iteration [15/20] [1987/1987]

Main iteration [16/20] [1987/1987]

Main iteration [17/20] [1987/1987]

Main iteration [18/20] [1987/1987]

Main iteration [19/20] [1987/1987]

Main iteration [20/20] [1987/1987]

Normalising graphs [1987/1987]

Solving haplotypes [1987/1987]

Writing sample list in [dog10k_plink_chr2.phased.sample]

Writing site list and haplotypes in [dog10k_plink_chr2.phased.haps]

Running time: 70320 seconds

Segmented HAPlotype Estimation & Imputation Tool
  * Authors : Olivier Delaneau, Jared O'Connell, Jean-François Zagury, Jonathan Marchini
  * Contact : send an email to the OXSTATGEN mail list https://www.jiscmail.ac.uk/cgi-bin/webadmin?A0=OXSTATGEN
  * Webpage : https://mathgen.stats.ox.ac.uk/shapeit
  * Version : v2.r904
  * Date    : 13/06/2025 07:43:08
  * LOGfile : [shapeit_13062025_07h43m08s_758ebcfc-65c3-40fe-9b7e-0202fbb5a15b.log]

MODE -convert : CONVERTING/SUBSETTING SHAPEIT HAPLOTYPES

Parameters :
  * Seed : 1749796988
  * Parallelisation: 1 threads
  * Ref allele is NOT aligned on the reference genome

Reading sample list in [dog10k_plink_chr2.phased.sample]
  * 1987 individuals included

Reading site list and haplotypes in [dog10k_plink_chr2.phased.haps]
  * 483737 sites included

Mapping sites and samples of the haplotype set
  * 483737 sites found
  * 0 sites not found
  * 3974 founder haplotypes

Writing sample list in [dog10k_plink_chr2.phased.impute.samples]

Writing site list in [dog10k_plink_chr2.phased.impute.legend]

Writing haplotypes in [dog10k_plink_chr2.phased.impute.haplotypes]

Running time: 474 seconds
```

</details>

---

### 5. Imputation

For each downsampled, prepared and pre-phased *Ostrander* chromosome, we want to now impute it against our new prepared and phased Dog10K reference panel.
For each Ostrander Downsampled Chromosome XX we now perform the imputation analysis using *impute2* against our *dog10k* reference.
Again we use a HPC script (below) and we also need a genetic map.
The [genetic map used](https://github.com/cflerin/dog_recombination)[^4] has a map of each Chromosome, this will be used for all phasing and imputation.

The imputation example below is carried out across the *entire* chromosome in each case we provide the chromosome length ($CHRSIZE) from the [karyotype file](etc/Canis_lupus_familiaris.ROS_Cfam_1.0.114.karyotype.tsv) described above.
  * The maps used are in the [maps folder](./maps).
  * ```impute2 -use_prephased_g -m ../maps/chr2.cf3.1_map.txt -h dog10k_plink_chr2.phased.impute.haplotypes -l dog10k_plink_chrXX.phased.impute.legend -known_haps_g broad_plink_chrXX_gwas.phased.haps -int 1 $CHRSIZE -allow_large_regions -Ne 200 -o broad_plink_chrXX_gwas.phased.impute_final -phase```

We can also break each chromosome into chunks of approximately 5Mb to 10Mb and run these individually, this has some key advantages:
  * It's supposedly more accurate as *Impute2* works better over regions of this size.
  * It's significantly faster as large chromosomes take a number of hours as an entire unit.
  * Later we will explore which option yields better results.

* Performing the imputation in 5Mb chunks would look like follows (first chunk):
  * `impute2 -use_prephased_g -m ../maps/chr2.cf3.1_map.txt -h dog10k_plink_chr2.phased.impute.haplotypes -l dog10k_plink_chrXX.phased.impute.legend -known_haps_g broad_plink_chrXX_gwas.phased.haps -int 1 5000000 -allow_large_regions -Ne 200 -o broad_plink_chrXX_gwas.phased.impute_final -phase`
  * the output results (per chunk) can then be concatenated.
  * The script below will use 5Mb chunking and create sub_folders for the results of each imputation, e.g. `impute_chr2/`.

### HPC Details:
* **Runtime - Approximately 10-28 hours per chromosome on 1 cpu each (whole chromsomomes).**
* **Runtime - Approximately up to 30 mins per chromosome run in parallel (5Mb chunks).**
* One script is used for this:
  * [impute.pl](scripts/impute.pl) - A perl script that generates an *sbatch* script for each chromosome that will launch via *sbatch* on our HPC.
* To modify the settings around whole chromosome or chunks edit the `$chunksize` option:
  * `$chunksize=5000000;` - Split each chromosome into 5Mb chunks
  * `$chunksize=1000000000000000000000;` - Set arbitrarily large value to force one chunk per chromosome.
  * The chunksize option is intelligent, it will work out an appropriate and equally sized number of chunks of a size close to, but slightly smaller, than that requested.
  * We may want to explore `-k` options and running full *MCMC* to improve accuracy later.

* run without parameters it will print the commands and scripts without submitting.
* To actually submit the jobs run with the `--runnit` flag, e.g. `./shape_it_gwas.pl --runnit`
  * sbatch parameters used:
    * `--partition=PlanEx2` (use the CGS cluster).
    * `--mem=60G` (60Gb RAM needed) - This can be much less (e.g. 3Gb for 5Mb chunks).
    * `--ntasks=1` (1 cpu per job)
    * `--cpus-per-task=1` (1 threads per job)
    * `-o XXXX,out` (output log)
    * `-e YYYY.err` (error log)

<details>
  <summary>Files produced for one Chromosome (chr2, 18 chunks of approx 5Mb)</summary>

```
broad_plink_chr2_gwas.phased.impute_final_chunk1
broad_plink_chr2_gwas.phased.impute_final_chunk10
broad_plink_chr2_gwas.phased.impute_final_chunk10_allele_probs
broad_plink_chr2_gwas.phased.impute_final_chunk10_haps
broad_plink_chr2_gwas.phased.impute_final_chunk10_info
broad_plink_chr2_gwas.phased.impute_final_chunk10_info_by_sample
broad_plink_chr2_gwas.phased.impute_final_chunk10_summary
broad_plink_chr2_gwas.phased.impute_final_chunk10_warnings
broad_plink_chr2_gwas.phased.impute_final_chunk11
broad_plink_chr2_gwas.phased.impute_final_chunk11_allele_probs
broad_plink_chr2_gwas.phased.impute_final_chunk11_haps
broad_plink_chr2_gwas.phased.impute_final_chunk11_info
broad_plink_chr2_gwas.phased.impute_final_chunk11_info_by_sample
broad_plink_chr2_gwas.phased.impute_final_chunk11_summary
broad_plink_chr2_gwas.phased.impute_final_chunk11_warnings
broad_plink_chr2_gwas.phased.impute_final_chunk12
broad_plink_chr2_gwas.phased.impute_final_chunk12_allele_probs
broad_plink_chr2_gwas.phased.impute_final_chunk12_haps
broad_plink_chr2_gwas.phased.impute_final_chunk12_info
broad_plink_chr2_gwas.phased.impute_final_chunk12_info_by_sample
broad_plink_chr2_gwas.phased.impute_final_chunk12_summary
broad_plink_chr2_gwas.phased.impute_final_chunk12_warnings
broad_plink_chr2_gwas.phased.impute_final_chunk13
broad_plink_chr2_gwas.phased.impute_final_chunk13_allele_probs
broad_plink_chr2_gwas.phased.impute_final_chunk13_haps
broad_plink_chr2_gwas.phased.impute_final_chunk13_info
broad_plink_chr2_gwas.phased.impute_final_chunk13_info_by_sample
broad_plink_chr2_gwas.phased.impute_final_chunk13_summary
broad_plink_chr2_gwas.phased.impute_final_chunk13_warnings
broad_plink_chr2_gwas.phased.impute_final_chunk14
broad_plink_chr2_gwas.phased.impute_final_chunk14_allele_probs
broad_plink_chr2_gwas.phased.impute_final_chunk14_haps
broad_plink_chr2_gwas.phased.impute_final_chunk14_info
broad_plink_chr2_gwas.phased.impute_final_chunk14_info_by_sample
broad_plink_chr2_gwas.phased.impute_final_chunk14_summary
broad_plink_chr2_gwas.phased.impute_final_chunk14_warnings
broad_plink_chr2_gwas.phased.impute_final_chunk15
broad_plink_chr2_gwas.phased.impute_final_chunk15_allele_probs
broad_plink_chr2_gwas.phased.impute_final_chunk15_haps
broad_plink_chr2_gwas.phased.impute_final_chunk15_info
broad_plink_chr2_gwas.phased.impute_final_chunk15_info_by_sample
broad_plink_chr2_gwas.phased.impute_final_chunk15_summary
broad_plink_chr2_gwas.phased.impute_final_chunk15_warnings
broad_plink_chr2_gwas.phased.impute_final_chunk16
broad_plink_chr2_gwas.phased.impute_final_chunk16_allele_probs
broad_plink_chr2_gwas.phased.impute_final_chunk16_haps
broad_plink_chr2_gwas.phased.impute_final_chunk16_info
broad_plink_chr2_gwas.phased.impute_final_chunk16_info_by_sample
broad_plink_chr2_gwas.phased.impute_final_chunk16_summary
broad_plink_chr2_gwas.phased.impute_final_chunk16_warnings
broad_plink_chr2_gwas.phased.impute_final_chunk17_summary
broad_plink_chr2_gwas.phased.impute_final_chunk17_warnings
broad_plink_chr2_gwas.phased.impute_final_chunk18
broad_plink_chr2_gwas.phased.impute_final_chunk18_allele_probs
broad_plink_chr2_gwas.phased.impute_final_chunk18_haps
broad_plink_chr2_gwas.phased.impute_final_chunk18_info
broad_plink_chr2_gwas.phased.impute_final_chunk18_info_by_sample
broad_plink_chr2_gwas.phased.impute_final_chunk18_summary
broad_plink_chr2_gwas.phased.impute_final_chunk18_warnings
broad_plink_chr2_gwas.phased.impute_final_chunk1_allele_probs
broad_plink_chr2_gwas.phased.impute_final_chunk1_haps
broad_plink_chr2_gwas.phased.impute_final_chunk1_info
broad_plink_chr2_gwas.phased.impute_final_chunk1_info_by_sample
broad_plink_chr2_gwas.phased.impute_final_chunk1_summary
broad_plink_chr2_gwas.phased.impute_final_chunk1_warnings
broad_plink_chr2_gwas.phased.impute_final_chunk2
broad_plink_chr2_gwas.phased.impute_final_chunk2_allele_probs
broad_plink_chr2_gwas.phased.impute_final_chunk2_haps
broad_plink_chr2_gwas.phased.impute_final_chunk2_info
broad_plink_chr2_gwas.phased.impute_final_chunk2_info_by_sample
broad_plink_chr2_gwas.phased.impute_final_chunk2_summary
broad_plink_chr2_gwas.phased.impute_final_chunk2_warnings
broad_plink_chr2_gwas.phased.impute_final_chunk3
broad_plink_chr2_gwas.phased.impute_final_chunk3_allele_probs
broad_plink_chr2_gwas.phased.impute_final_chunk3_haps
broad_plink_chr2_gwas.phased.impute_final_chunk3_info
broad_plink_chr2_gwas.phased.impute_final_chunk3_info_by_sample
broad_plink_chr2_gwas.phased.impute_final_chunk3_summary
broad_plink_chr2_gwas.phased.impute_final_chunk3_warnings
broad_plink_chr2_gwas.phased.impute_final_chunk4
broad_plink_chr2_gwas.phased.impute_final_chunk4_allele_probs
broad_plink_chr2_gwas.phased.impute_final_chunk4_haps
broad_plink_chr2_gwas.phased.impute_final_chunk4_info
broad_plink_chr2_gwas.phased.impute_final_chunk4_info_by_sample
broad_plink_chr2_gwas.phased.impute_final_chunk4_summary
broad_plink_chr2_gwas.phased.impute_final_chunk4_warnings
broad_plink_chr2_gwas.phased.impute_final_chunk5
broad_plink_chr2_gwas.phased.impute_final_chunk5_allele_probs
broad_plink_chr2_gwas.phased.impute_final_chunk5_haps
broad_plink_chr2_gwas.phased.impute_final_chunk5_info
broad_plink_chr2_gwas.phased.impute_final_chunk5_info_by_sample
broad_plink_chr2_gwas.phased.impute_final_chunk5_summary
broad_plink_chr2_gwas.phased.impute_final_chunk5_warnings
broad_plink_chr2_gwas.phased.impute_final_chunk6
broad_plink_chr2_gwas.phased.impute_final_chunk6_allele_probs
broad_plink_chr2_gwas.phased.impute_final_chunk6_haps
broad_plink_chr2_gwas.phased.impute_final_chunk6_info
broad_plink_chr2_gwas.phased.impute_final_chunk6_info_by_sample
broad_plink_chr2_gwas.phased.impute_final_chunk6_summary
broad_plink_chr2_gwas.phased.impute_final_chunk6_warnings
broad_plink_chr2_gwas.phased.impute_final_chunk7
broad_plink_chr2_gwas.phased.impute_final_chunk7_allele_probs
broad_plink_chr2_gwas.phased.impute_final_chunk7_haps
broad_plink_chr2_gwas.phased.impute_final_chunk7_info
broad_plink_chr2_gwas.phased.impute_final_chunk7_info_by_sample
broad_plink_chr2_gwas.phased.impute_final_chunk7_summary
broad_plink_chr2_gwas.phased.impute_final_chunk7_warnings
broad_plink_chr2_gwas.phased.impute_final_chunk8
broad_plink_chr2_gwas.phased.impute_final_chunk8_allele_probs
broad_plink_chr2_gwas.phased.impute_final_chunk8_haps
broad_plink_chr2_gwas.phased.impute_final_chunk8_info
broad_plink_chr2_gwas.phased.impute_final_chunk8_info_by_sample
broad_plink_chr2_gwas.phased.impute_final_chunk8_summary
broad_plink_chr2_gwas.phased.impute_final_chunk8_warnings
broad_plink_chr2_gwas.phased.impute_final_chunk9
broad_plink_chr2_gwas.phased.impute_final_chunk9_allele_probs
broad_plink_chr2_gwas.phased.impute_final_chunk9_haps
broad_plink_chr2_gwas.phased.impute_final_chunk9_info
broad_plink_chr2_gwas.phased.impute_final_chunk9_info_by_sample
broad_plink_chr2_gwas.phased.impute_final_chunk9_summary
broad_plink_chr2_gwas.phased.impute_final_chunk9_warnings
```

</details>



<details>
  <summary>Output Imputation Log Example (Chromosome 2)</summary>

```

======================
 IMPUTE version 2.3.2 
======================

Copyright 2008 Bryan Howie, Peter Donnelly, and Jonathan Marchini
Please see the LICENCE file included with this program for conditions of use.

The seed for the random number generator is 1632239422.

Command-line input: /mnt/research2/anton/Eleanor_Raffan/impute_v2.3.2_x86_64_static/impute2 -use_prephased_g -m /mnt/research2/anton/Eleanor_Raffan/maps/chr2.cf3.1_map.txt -h /mnt/research2/anton/Eleanor_Raffan/dog10k/panel/dog10k_plink_chr2.phased.impute.haplotypes -l /mnt/research2/anton/Eleanor_Raffan/dog10k/panel/dog10k_plink_chr2.phased.impute.legend -known_haps_g broad_plink_chr2_gwas.phased.haps -int 1 4788212 -allow_large_regions -Ne 200 -o broad_plink_chr2_gwas.phased.impute_final_chunk1 -phase

---------------------------------
 Nomenclature and data structure 
---------------------------------

     Panel 0: phased reference haplotypes
     Panel 2: phased study haplotypes

For optimal results, each successive panel (0,1,2) should contain a subset of the SNPs in the previous panel. When the data structure deviates from this ideal configuration, IMPUTE2 tries to use as much of the available information as possible; see documentation for details.

-------------
 Input files 
-------------

         Panel 0 haplotypes: /mnt/research2/anton/Eleanor_Raffan/dog10k/panel/dog10k_plink_chr2.phased.impute.haplotypes
         Panel 0 hap legend: /mnt/research2/anton/Eleanor_Raffan/dog10k/panel/dog10k_plink_chr2.phased.impute.legend
         Panel 2 known haps: broad_plink_chr2_gwas.phased.haps
                genetic map: /mnt/research2/anton/Eleanor_Raffan/maps/chr2.cf3.1_map.txt

--------------
 Output files 
--------------

                main output: broad_plink_chr2_gwas.phased.impute_final_chunk1
                SNP QC info: broad_plink_chr2_gwas.phased.impute_final_chunk1_info
             sample QC info: broad_plink_chr2_gwas.phased.impute_final_chunk1_info_by_sample
                run summary: broad_plink_chr2_gwas.phased.impute_final_chunk1_summary
                warning log: broad_plink_chr2_gwas.phased.impute_final_chunk1_warnings
        Panel 2 phased haps: broad_plink_chr2_gwas.phased.impute_final_chunk1_haps
       Panel 2 allele probs: broad_plink_chr2_gwas.phased.impute_final_chunk1_allele_probs

-----------------
 Data processing 
-----------------

-reading genetic map from -m file
 --filename=[/mnt/research2/anton/Eleanor_Raffan/maps/chr2.cf3.1_map.txt]
 --read 5 SNPs in the analysis interval+buffer region

-reading Panel 2 haplotypes from -known_haps_g file
 --filename=[broad_plink_chr2_gwas.phased.haps]
 --detected 668 individuals
 --read 116 SNPs in the analysis interval+buffer region
 --added 116 new SNPs based on known haplotypes

-reading Panel 0 haplotypes from -h and -l files
 --filename=[/mnt/research2/anton/Eleanor_Raffan/dog10k/panel/dog10k_plink_chr2.phased.impute.haplotypes]
 --filename=[/mnt/research2/anton/Eleanor_Raffan/dog10k/panel/dog10k_plink_chr2.phased.impute.legend]
 --detected 3974 haplotypes
 --read 17889 SNPs in the analysis interval+buffer region

-removing SNPs that violate the hierarchical data requirements
 --no SNPs removed

-removing reference-only SNPs from buffer region
 --removed 791 SNPs

-checking strand alignment between Panel 2 and Panel 0 by allele labels
 --flipped strand due to allele mismatch at 0 out of 116 SNPs in Panel 2

-aligning allele labels between panels

-removing non-aligned genotyped SNPs
 --removed 0 out of 112 SNPs with data in multiple panels

--------------
 Data summary 
--------------

[type 0 = SNP in Panel 0 only]
[type 1 = SNP in Panel 1]
[type 2 = SNP in Panel 2 and all ref panels]
[type 3 = SNP in Panel 2 only]

-Upstream buffer region
 --0 type 0 SNPs
 --0 type 1 SNPs
 --0 type 2 SNPs
 --0 type 3 SNPs
 --0 total SNPs

-Downstream buffer region
 --0 type 0 SNPs
 --0 type 1 SNPs
 --7 type 2 SNPs
 --0 type 3 SNPs
 --7 total SNPs

-Analysis region (as defined by -int argument)
 --16986 type 0 SNPs
 --0 type 1 SNPs
 --105 type 2 SNPs
 --4 type 3 SNPs
 --17095 total SNPs

-Output file
 --16986 type 0 SNPs
 --0 type 1 SNPs
 --105 type 2 SNPs
 --4 type 3 SNPs

-In total, 17102 SNPs will be used in the analysis, including 112 Panel 2 SNPs

-setting storage space

----------------
 Run parameters 
----------------

        reference haplotypes: 3974 [Panel 0]
           study individuals: 668 [Panel 2]
           sequence interval: [1,4788212]
                      buffer: 250 kb
                          Ne: 200
           input call thresh: 0.900
     burn-in MCMC iterations: 0
       total MCMC iterations: 1 (1 used for inference)
   HMM states for imputation: 500 [Panel 0->2]
                active flags: <-use_prephased_g> <-allow_large_regions> <-phase>

---------
 Run log 
---------

RESETTING PARAMETERS FOR "SURROGATE FAMILY" MODELING
-setting mutation matrices
-setting switch rates

MCMC iteration [1/1]


diploid sampling success rate: (no diploid sampling performed)

haploid sampling success rate: (no haploid sampling performed)


--------------------------------
 Imputation accuracy assessment 
--------------------------------

The table below is based on an internal cross-validation that is performed during each IMPUTE2 run. For this analysis, the program masks the genotypes of one variant at a time in the study data (Panel 2) and imputes the masked genotypes by using the remaining study and reference data. The imputed genotypes are then compared with the original genotypes to produce the concordance statistics shown in the table. You can learn more about this procedure and the contents of the table at http://mathgen.stats.ox.ac.uk/impute/concordance_table_description.html.

In the current analysis, IMPUTE2 masked, imputed, and evaluated 70140 genotypes that were called with high confidence (maximum probability >= 0.90) in the Panel 2 input file (-g or -known_haps_g).

When the masked study genotypes were imputed with reference data from Panel 0, the concordance between original and imputed genotypes was as follows:

  Interval  #Genotypes %Concordance         Interval  %Called %Concordance
  [0.0-0.1]          0          0.0         [ >= 0.0]   100.0         97.3
  [0.1-0.2]          0          0.0         [ >= 0.1]   100.0         97.3
  [0.2-0.3]          0          0.0         [ >= 0.2]   100.0         97.3
  [0.3-0.4]          0          0.0         [ >= 0.3]   100.0         97.3
  [0.4-0.5]        141         30.5         [ >= 0.4]   100.0         97.3
  [0.5-0.6]        649         59.0         [ >= 0.5]    99.8         97.4
  [0.6-0.7]        643         68.4         [ >= 0.6]    98.9         97.8
  [0.7-0.8]        843         74.1         [ >= 0.7]    98.0         98.1
  [0.8-0.9]       1264         84.7         [ >= 0.8]    96.8         98.4
  [0.9-1.0]      66600         98.6         [ >= 0.9]    95.0         98.6

-generating consensus haplotype estimates (relative to fixed input haps)

Have a nice day!
```
</details>


---

# Imputation Concordance Analysis
We now use the imputed dataset to explore QC metrics for imputation, for variants globally, per chromosome and  across individuals.
We can do this two ways, buy using the Impute2 QC metrics and also comparing back to the original truth dataset from WGS.


#### Overview of Concordance Process

For the two ways to assess this:

#### 1. Exploration of Impute2 Statistics
* **Impute2** allows robust exploration of imputation quality.
* For SNPs we supplied from Ostrander, the input genotypes at that SNP were  masked internally and then imputed as if the SNP were of type X. Similarly, r2_typeX is the squared correlation between input and masked/imputed genotypes at a SNP. This allows quite a robust way for Impute 2 to assess imputation accuracy.
* We have 2 values we can explore:
  * Concordance: *concord_type0* 
  * $r^2$: - *r2_type0*
* We can explore these by SNP Position, Chromosome Chunk, Chromosome.
* Additionally, we also have sample specific info files so can explore how different individual dogs/breeds/types fare.

#### Example QC output from *impute2*
This is the type of data produced for chromosome1. The last two columns show the concordance and $r^2$.
```
snp_id rs_id position a0 a1 exp_freq_a1 info certainty type info_type0 concord_type0 r2_type0
1 1:34533898 34533898 T C 0.824 1.000 1.000 2 0.876 0.896 0.714
1 1:34556174 34556174 C A 0.497 1.000 1.000 2 0.923 0.885 0.786
1 1:34582419 34582419 A C 0.584 1.000 1.000 2 0.922 0.896 0.808
1 1:34598722 34598722 C T 0.584 1.000 1.000 2 0.906 0.867 0.745
1 1:34702504 34702504 A G 0.665 1.000 1.000 2 0.883 0.826 0.711
1 1:34730426 34730426 T C 0.715 1.000 1.000 2 0.883 0.836 0.648
1 1:34742583 34742583 C T 0.762 1.000 1.000 2 0.903 0.863 0.678
1 1:34766687 34766687 T C 0.546 1.000 1.000 2 0.904 0.846 0.747
1 1:34782808 34782808 C T 0.691 1.000 1.000 2 0.878 0.851 0.713
1 1:34800994 34800994 T C 0.796 1.000 1.000 2 0.864 0.845 0.620
1 1:34831956 34831956 A G 0.806 1.000 1.000 2 0.880 0.875 0.649
1 1:34883307 34883307 A G 0.592 1.000 1.000 2 0.888 0.875 0.778
1 1:34908754 34908754 A G 0.817 1.000 1.000 2 0.978 0.984 0.960
1 1:34919792 34919792 T G 0.835 1.000 1.000 2 0.981 0.988 0.974
1 1:34946032 34946032 A G 0.805 1.000 1.000 2 0.977 0.976 0.935
1 1:34956621 34956621 C A 0.702 1.000 1.000 2 0.944 0.934 0.882
1 1:35008786 35008786 G T 0.746 1.000 1.000 2 0.928 0.918 0.825
1 1:35035399 35035399 A G 0.921 1.000 1.000 2 0.947 0.964 0.828
```

We also have QC metrics (concordance and $r^2$ for each individual in that chromosome chunk.
In this case each individual has a specific line in their original order.
```
concord_type0 r2_type0
0.986 0.977
0.966 0.943
1.000 1.000
0.979 0.950
0.973 0.939
0.986 0.970
0.945 0.915
0.945 0.889
0.884 0.751
0.842 0.557
0.973 0.946
0.986 0.977
1.000 0.996
```

#### Exploration of *Impute2* QC metrics.

We use a number of perl scripts to globally explore the concordance and $r^2$ values:
* Median, Standard Deviation, Average, Max and Min for:
* Whole Chromosome chunks
* Whole Chromosomes
* Individuals by Chromosome
* All Individuals globally

These scripts generate *txt* files we can explore in *R/BioConductor*.

![histogram]("analysis/Imputation_Analysis_files/figure-gfm/unnamed-chunk-2-1.png")


#### 2. Evaluation of raw SNP non-reference concordance back to truth dataset.

To perform the concordance analysis we work as follows using the [concordance.pl](scripts/concordance.pl) script:

#### Concordance Analysis Script
Invocation is as follows for a single chromosome (XX):

`./concordance.pl broad_plink_chrXX_gwas.phased.impute_final_haps.merged ../../broad-chr2.cf4.vcf.gz`

This will create two new files:
* `broad_plink_chrXX_gwas.phased.impute_final_haps.merged.concordance_ind.txt`
* `broad_plink_chrXX_gwas.phased.impute_final_haps.merged.concordance_pos.txt`
* Process the imputation haplotype file for each chromosome and store the REF and ALT allele and haplotype calls.
* Process the imputation haplotype file and store the MAF for each position.
* Extract the sample list of 676 dogs that passed filtering through the dog10k process.
* Extract the sample list of 668 dogs from the Ostrander VCF files
* Process the Ostrander chromosome (liftover) VCF file, SNV by SNV
  * If a position matches:
    * iterate over individuals that match and call their haplotypes on both and see if they match
    * Store concordant counts and disconcordant counts for all assessible positions.
    * Compute the percentage of comparable positions that are concordant
    * We could also compute an $r^2$.
  * Save the concordance data globally, per dog and per snp to a new file.
* We would now be able to explore the concordance in multiple ways.
  * Global concordance.
  * Concordance vs MAF.
  * Concordance by Dog/Breed.
  * Concordance by position on Chromosome.

---

## References
[^1]:Tong Zhou, Shao-Yan Pu, Shao-Jie Zhang, Qi-Jun Zhou, Min Zeng, Jing-Sheng Lu, Xuemei Lu, Ya-Nan Wang, Guo-Dong Wang, Dog10K: an integrated Dog10K database summarizing canine multi-omics, Nucleic Acids Research, V53(D1), 6 Jan (2025).
[^2]:Plassais J, Kim J, Davis BW, Karyadi DM, Hogan AN, Harris AC, Decker B, Parker HG, Ostrander EA. Whole genome sequencing of canids reveals genomic regions under selection and variants influencing morphology. Nat Commun. 10:1489 (2019).
[^3]:https://www.wisdompanel.com/en-gb.
[^4]:Campbell, CL, Bhérer, C, Morrow, BE, Boyko, AR, & Auton, A. A Pedigree-Based Map of Recombination in the Domestic Dog Genome. G3: Genes, Genomes, Genetics, v6 no11 3517-3524; (2016)
[^5]:Genovese G., McCarroll S. et al. BCFtools/liftover: an accurate and comprehensive tool to convert genetic variants across genome assemblies. Bioinformatics 40;2 (2024).

