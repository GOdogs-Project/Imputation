#!/bin/sh
echo $1
echo $2

plink --const-fid 0 --vcf ../broad-chr$2.downsampled.cf4.vcf.gz --out broad_plink_chr$2 --dog
perl -lane 'next if($F[0] eq "0"); print $F[1]."\t".$F[0].":".$F[3];' broad_plink_chr$2.bim > broad_plink_chr$2.names
plink --update-name broad_plink_chr$2.names --make-bed --bfile broad_plink_chr$2 --out broad_plink_chr$2.1 --dog
perl -lane 'if(($F[4] eq "A" && $F[5] eq "T") || ($F[4] eq "T" && $F[5] eq "A") || ($F[4] eq "C" && $F[5] eq "G") || ($F[4] eq "G" && $F[5] eq "C")){ print $F[1] }' broad_plink_chr$2.1.bim > broad_plink_ambiguous.snps.chr$2
plink --bfile broad_plink_chr$2.1 --exclude broad_plink_ambiguous.snps.chr$2 --make-bed --out broad_plink_chr$2.2 --dog
plink --bfile broad_plink_chr$2.2 --geno 0.03 --make-bed --out broad_plink_chr$2.geno --dog
plink --bfile broad_plink_chr$2.geno --mind 0.1 --make-bed --out broad_plink_chr$2.mind --dog
plink --bfile broad_plink_chr$2.mind --maf 0.01 --make-bed --out broad_plink_chr$2.maf --dog
plink --bfile broad_plink_chr$2.maf --hwe 0.00005 --make-bed --out broad_plink_chr$2_gwas --dog
