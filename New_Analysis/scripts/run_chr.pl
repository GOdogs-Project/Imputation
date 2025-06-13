#!/usr/bin/perl

$batch="sbatch --partition=PlanEx2 --mem=2G --ntasks 1 -o XXXX -e YYYY";

foreach $file (@ARGV){
	chomp($file);
	$stem=$file;
	$stem=~ s/^.*chr(\d+).*$/$1/g;
	$cmd="srun --mem=2G ./run_chr.sh $file $stem";
	open(FILEOUT,"> /tmp/sbatch$$\_chr$stem.sh");
	print FILEOUT "#!/bin/sh\n$cmd\n";
	close(FILEOUT);

	$sbatch=$batch;
	$sbatch=~ s/XXXX/chr$stem.out/g;
	$sbatch=~ s/YYYY/chr$stem.err/g;

	print "$sbatch /tmp/sbatch$$\_chr$stem.sh\n";
	system("$sbatch /tmp/sbatch$$\_chr$stem.sh");
}

