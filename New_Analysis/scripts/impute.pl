#!/usr/bin/perl
use POSIX;

$impute="/mnt/research2/anton/Eleanor_Raffan/impute_v2.3.2_x86_64_static/impute2";
$maps="/mnt/research2/anton/Eleanor_Raffan/maps";
$dog10k="/mnt/research2/anton/Eleanor_Raffan/dog10k/panel";
$rcmd="$impute -use_prephased_g -m $maps/chrXX.cf3.1_map.txt -h $dog10k/dog10k_plink_chrXX.phased.impute.haplotypes -l $dog10k/dog10k_plink_chrXX.phased.impute.legend -known_haps_g broad_plink_chrXX_gwas.phased.haps -int ZZ WW  -allow_large_regions -Ne 200 -o broad_plink_chrXX_gwas.phased.impute_final_chunkPP -phase";
$karyotype = "/mnt/research2/anton/Eleanor_Raffan/BroadInstitute_Ostrander_genomes/Canis_lupus_familiaris.ROS_Cfam_1.0.114.karyotype.tsv";

$sbatch="sbatch --partition=PlanEx2 --mem=10G --cpus-per-task 1 --ntasks 1 -o XX -e YY";
$chunksize=5000000;

# Setting this will effectively disable chromosome chunking
#$chunksize=1000000000000000000000;

$run=0;

foreach $arg(@ARGV){
	if ($arg eq '--runnit'){
		$run=1;
	}
}

open(FILE,$karyotype);
while(<FILE>){
	chomp;
	@array=split("\t",$_);
	$chr_size{$array[0]}=$array[1];

}

#for ($i=1;$i<=38;$i++){
for ($i=2;$i<=2;$i++){
	print "CHROMOSOME:$i CHRSIZE: $chr_size{$i}\n";
	$chunks=ceil($chr_size{$i}/$chunksize);
	print "Will divide into $chunks chunks\n";
	$size=ceil($chr_size{$i}/$chunks);

	print "Will divide into $chunks chunks of size $size bp\n";
	for ($j=0;$j<$chunks;$j++){
		$outfile="/tmp/sbatch_impute_$$\_$i\_chunk".($j+1).".sh";
		open(FILEOUT,">$outfile");
        	print FILEOUT "#!/bin/sh\n";
		$start=($j*$size);
		$end=$start+$size;
		if ($end > $chr_size{$i}){
			$end=$chr_size{$i};
		}
		$start=$start+1;
		$cmd=$rcmd;
		$this_chunk=$j+1;
        	$cmd=~ s/XX/$i/g;
		$cmd=~ s/ZZ/$start/g;
		$cmd=~ s/WW/$end/g;
		$cmd=~ s/PP/$this_chunk/g;
		print FILEOUT "$cmd\n";
		print "BATCH IN: $outfile\n";
		$batch=$sbatch;
		$batch=~ s/XX/sbatch_impute_chr$i\_chunk$this_chunk.out/g;
		$batch=~ s/YY/sbatch_impute_chr$i\_chunk$this_chunk.err/g;
		print "$batch $outfile\n";
		if ($run){
			system("$batch $outfile");
		}
		close FILEOUT;
		print "\n\n";
	}
}

