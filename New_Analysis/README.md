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
|Dog10K | 51M  | 1987 | CanFam 4            |
|Ostrander| 91M | 676 | CanFam 3.1          |
|Wisdom|95,165 |  | CanFam 3.1              |


## Software Used

* [plink](https://www.cog-genomics.org/plink/) - PLINK v1.90b6.24 64-bit (6 Jun 2021)
* [plink2](https://www.cog-genomics.org/plink/2.0/) - PLINK v2.0.0-a.6.16LM AVX2 Intel (9 Jun 2025)
* [bcftools](https://github.com/samtools/bcftools) - 1.19-46-gc63329bb
* [samtools](https://github.com/samtools/samtools) - 1.19.2-7-g0c03cbd
* [score](https://github.com/freeseek/score) - v1 2024[^5]
* [shapeit](https://mathgen.stats.ox.ac.uk/shapeit) - Version : v2.r904 
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
* `./downsample.pl broad-chr*.vcf.gz`
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
* Downsampling is run by a custom perl script [downsample.pl](scripts/downsample.pl).
  * Invocation: `./downsample.pl *downsampled.cf4.vcf.gz`

* Prepare a PLINK dataset for each downsampled Chromosome (XX where XX is Chr1 to Chr38).
  * `plink --const-fid 0 --vcf ../broad-chrXX.downsampled.cf4.vcf.gz --out broad_plink_chrXX --dog`
* Fix the Name column as it still contains CHR:POS names from *cf3.1*:
  * `perl -lane 'next if($F[0] eq "0"); print $F[1]."\t".$F[0].":".$F[3];' broad_plink_chrXX.bim > broad_plink_chrXX.names`
  * `plink --update-name broad_plink_chrXX.names --make-bed --bfile broad_plink_chrXX --out broad_plink_chrXX.1 --dog`
* Find and store ambiguous SNPs from the new plink file
  * `perl -lane 'if(($F[4] eq "A" && $F[5] eq "T") || ($F[4] eq "T" && $F[5] eq "A") || ($F[4] eq "C" && $F[5] eq "G") || ($F[4] eq "G" && $F[5] eq "C")){ print $F[1] }' broad_plink_chrXX.1.bim > ambiguous.snps.chrXX`
* Exclude the ambiguous SNPs in a new plink file
  * `plink --bfile broad_plink_chrXX.1 --exclude ambiguous.snps.chrXX --make-bed --out broad_plink_chrXX.2 --dog`
* Continue with filtering
  * `plink --bfile broad_plink_chrXX.2 --geno 0.03 --make-bed --out broad_plink_chrXX.geno --dog`
  * `plink --bfile broad_plink_chrXX.geno --mind 0.1 --make-bed --out broad_plink_chrXX.mind --dog`
  * `plink --bfile broad_plink_chrXX.mind --maf 0.01 --make-bed --out broad_plink_chrXX.maf --dog`
  * `plink --bfile broad_plink_chrXX.maf --hwe 0.00005 --make-bed --out broad_plink_chrXX_gwas --dog`

### HPC Details:
* **Runtime - Downsampling, 5-10 mins per chromosome, PLINK: approximately 30-60secs per chromosome**.
* Two scripts are used for this:
  * [run_chr.sh](scripts/run_chr.sh) - A generic script that runs a downsampled chromosome through the plink commands above.
  * [run_chr.pl](scripts/run_chr.pl) - A perl script that generates an *sbatch* script for each chromosome that will launch *run_chr.sh* via *sbatch* on our HPC.
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

---

### 2. Ostrander Phasing

For each Chromosome XX we now perform the phasing analysis using *shapeit*, again we use a HPC script (below) and we also need a genetic map.
The [genetic map used](https://github.com/cflerin/dog_recombination)[^4] has a map of each Chromosome, this will be used for all phasing and imputation.
  * The maps used are in the [maps folder](./maps).
  * `shapeit -M ../../maps/chrXX.cf3.1_map.txt -B broad_plink_chrXX_gwas -O broad_plink_chrXX_gwas.phased  -T 4 --window 2 --effective-size 200 --force`
  * `shapeit -convert --input-haps broad_plink_chrXX_gwas.phased --output-ref broad_plink_chrXX_gwas.phased.impute`

### HPC Details:
* **Runtime - Approximately 2 mins per chromosome.**
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


---

### 3. Dog10K Preparation

We will take the giant PLINK formatted datafile `dog10k.SNPs.plink` from the Dog10K dataset and split it first into individual chromosome (XX is the Chr number), PLINK files.
We will also do some simple filtering and create a new *refpanel* plink object for each.
  * `plink --bfile ../dog10k.SNPs.plink --chr XX --make-bed --out dog10k_plink_chrXX --dog`
  * `plink --bfile dog10k_plink_chrXX --maf 0.01 --mind 0.1 --geno 0.03 --make-bed --out dog10k_plink_chrXX.refpanel --dog`

### HPC Details:
* **Runtime - Approximately 10 mins per chromosome.**
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
* **Runtime - Approximately 10-28 hours per chromosome on 20 cpus each.**
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
The imputation here is carried out across the *entire* chromosome in each case we provide the chromosome length ($CHRSIZE) from the [karyotype file](etc/Canis_lupus_familiaris.ROS_Cfam_1.0.114.karyotype.tsv) described above.
  * The maps used are in the [maps folder](./maps).
  * ```impute2 -use_prephased_g -m ../maps/chr2.cf3.1_map.txt -h dog10k_plink_chr2.phased.impute.haplotypes -l dog10k_plink_chrXX.phased.impute.legend -known_haps_g broad_plink_chrXX_gwas.phased.haps -int 1 $CHRSIZE -allow_large_regions -Ne 200 -o broad_plink_chrXX_gwas.phased.impute_final -phase```

### HPC Details:
* **Runtime - Approximately 10-28 hours per chromosome on 20 cpus each.**
* One script is used for this:
  * [impute.pl](scripts/impute.pl) - A perl script that generates an *sbatch* script for each chromosome that will launch via *sbatch* on our HPC.
  * run without parameters it will print the commands and scripts without submitting.
  * To actually submit the jobs run with the `--runnit` flag, e.g. `./shape_it_gwas.pl --runnit`
  * sbatch parameters used:
    * `--partition=PlanEx2` (use the CGS cluster).
    * `--mem=60G` (60Gb RAM needed)
    * `--ntasks=1` (1 cpu per job)
    * `--cpus-per-task=1` (1 threads per job)
    * `-o XXXX,out` (output log)
    * `-e YYYY.err` (error log)

---

# Concordance Analysis
We now use the imputed dataset to calculate *Non-Reference Concordance* (NRC) for each variant in each individual based on the original truth dataset from WGS.

---

## References
[^1]:Tong Zhou, Shao-Yan Pu, Shao-Jie Zhang, Qi-Jun Zhou, Min Zeng, Jing-Sheng Lu, Xuemei Lu, Ya-Nan Wang, Guo-Dong Wang, Dog10K: an integrated Dog10K database summarizing canine multi-omics, Nucleic Acids Research, V53(D1), 6 Jan (2025).
[^2]:Plassais J, Kim J, Davis BW, Karyadi DM, Hogan AN, Harris AC, Decker B, Parker HG, Ostrander EA. Whole genome sequencing of canids reveals genomic regions under selection and variants influencing morphology. Nat Commun. 10:1489 (2019).
[^3]:https://www.wisdompanel.com/en-gb.
[^4]:Campbell, CL, Bhérer, C, Morrow, BE, Boyko, AR, & Auton, A. A Pedigree-Based Map of Recombination in the Domestic Dog Genome. G3: Genes, Genomes, Genetics, v6 no11 3517-3524; (2016)
[^5]:Genovese G., McCarroll S. et al. BCFtools/liftover: an accurate and comprehensive tool to convert genetic variants across genome assemblies. Bioinformatics 40;2 (2024).

