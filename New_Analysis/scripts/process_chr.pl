#!/usr/bin/perl
#
# proces_chr.pl - Anton Enright (aje39@cam.ac.uk).
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

$run=0;
$batch="sbatch --partition=PlanEx2 --mem=10G --ntasks 1 -o XXXX -e YYYY";
foreach $arg(@ARGV){
	if ($arg eq '--runnit'){
		$run=1;
	}
}

#$cmd="srun --mem=8G plink --bfile ../dog10k.SNPs.plink --chr XXXX --make-bed --out YYYY --dog";
$cmd="plink --bfile ../dog10k.SNPs.plink --chr XXXX --make-bed --out YYYY --dog";
$cmd2="plink --bfile YYYY --maf 0.01 --mind 0.1 --geno 0.03 --make-bed --out YYYY.refpanel --dog; ";
$cmd3="/mnt/research2/anton/Eleanor_Raffan/plink2 --bfile YYYY.refpanel --set-all-var-ids \@:\# --make-bed --dog --out YYYY.refpanel.names;\n/mnt/research2/anton/Eleanor_Raffan/plink2 --bfile YYYY.refpanel.names --freq --dog --out YYYY.refpanel.names.maf";

for ($i=1;$i<=38;$i++){

	$rcmd=$cmd;
	$rcmd=~ s/XXXX/$i/g;
	$rcmd=~ s/YYYY/dog10k_plink_chr$i/g;
	$rcmd2=$cmd2;
	$rcmd2=~s/YYYY/dog10k_plink_chr$i/g;
	$rcmd3=$cmd3;
	$rcmd3=~s/YYYY/dog10k_plink_chr$i/g;
	print "$rcmd\n";
	open(FILE,">/tmp/$$\_sbatch$i.sh");
        print FILE "#!/bin/sh\n$rcmd\n$rcmd2\n$rcmd3\n";
	$sbatch=$batch;
	$sbatch=~ s/XXXX/sbatch_$i.out/g;
	$sbatch=~ s/YYYY/sbatch_$i.err/g;
	print "$sbatch /tmp/$$\_sbatch$i.sh\n";
	if ($run){
		system("$sbatch /tmp/$$\_sbatch$i.sh");
	}
}
