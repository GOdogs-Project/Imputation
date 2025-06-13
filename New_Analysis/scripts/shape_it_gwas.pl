#!/usr/bin/perl
$run=0;
$batch="sbatch --partition=PlanEx2 --mem=30G --cpus-per-task 4 --ntasks 1 -o XX -e YY";
foreach $arg(@ARGV){
	if ($arg eq '--runnit'){
		$run=1;
	}
}

$cmd="../../shapeit.v2.904.3.10.0-693.11.6.el7.x86_64/bin/shapeit -M ../../maps/chrXX.cf3.1_map.txt -B broad_plink_chrXX_gwas -O broad_plink_chrXX_gwas.phased -T 4 --window 2 --effective-size 200 --force\n../../shapeit.v2.904.3.10.0-693.11.6.el7.x86_64/bin/shapeit -convert --input-haps broad_plink_chrXX_gwas.phased --output-ref broad_plink_chrXX_gwas.phased.impute";
for ($i=1;$i<=38;$i++){

	$rcmd=$cmd;
	$rcmd=~ s/XX/$i/g;
	print "$rcmd\n";
	open(FILE,">/tmp/$$\_sbatch$i.sh");
        print FILE "#!/bin/sh\n$rcmd\n";
	$sbatch=$batch;
	$sbatch=~ s/XX/sbatch_shapeit_$i.out/g;
	$sbatch=~ s/YY/sbatch_shapeit_$i.err/g;
	print "$sbatch /tmp/$$\_sbatch$i.sh\n";
	if ($run){
		system("$sbatch /tmp/$$\_sbatch$i.sh");
	}
}

