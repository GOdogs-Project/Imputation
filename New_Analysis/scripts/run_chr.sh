#!/bin/sh
#
# run_chr.sh - Anton Enright (aje39@cam.ac.uk).
#
# Copyright (c) 2025 - Anton Enright - University of Cambridge
#
# This program is free software: you can redistribute it and/or modify it under the terms of the GNU General Public License as
# published by the Free Software Foundation, either version 3 of the License, or (at your option) any later version.
# This program is distributed in the hope that it will be useful, but WITHOUT ANY WARRANTY; without even the implied warranty 
# of MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE. See the GNU General Public License for more details.
#
# You should have received a copy of the GNU General Public License along with this program. If not, see <https://www.gnu.org/licenses/>.
#
echo $1
echo $2

plink --const-fid 0 --vcf $1 --chr $2 --out broad_plink_chr$2 --allow-extra-chr --dog
perl -lane 'next if($F[0] eq "0"); print $F[1]."\t".$F[0].":".$F[3];' broad_plink_chr$2.bim > broad_plink_chr$2.names
plink --update-name broad_plink_chr$2.names --make-bed --bfile broad_plink_chr$2 --out broad_plink_chr$2.1 --dog
perl -lane 'if(($F[4] eq "A" && $F[5] eq "T") || ($F[4] eq "T" && $F[5] eq "A") || ($F[4] eq "C" && $F[5] eq "G") || ($F[4] eq "G" && $F[5] eq "C")){ print $F[1] }' broad_plink_chr$2.1.bim > broad_plink_ambiguous.snps.chr$2
plink --bfile broad_plink_chr$2.1 --exclude broad_plink_ambiguous.snps.chr$2 --make-bed --out broad_plink_chr$2.2 --dog
plink --bfile broad_plink_chr$2.2 --geno 0.03 --make-bed --out broad_plink_chr$2.geno --dog
plink --bfile broad_plink_chr$2.geno --mind 0.1 --make-bed --out broad_plink_chr$2.mind --dog
plink --bfile broad_plink_chr$2.mind --maf 0.01 --make-bed --out broad_plink_chr$2.maf --dog
plink --bfile broad_plink_chr$2.maf --hwe 0.000000000000000000005 --make-bed --out broad_plink_chr$2_gwas --dog
