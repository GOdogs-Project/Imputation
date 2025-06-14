#!/usr/bin/perl
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
