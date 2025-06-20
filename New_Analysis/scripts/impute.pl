#!/usr/bin/perl
#
# impute.pl - Anton Enright (aje39@cam.ac.uk).
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
use POSIX;

$impute="/mnt/research2/anton/Eleanor_Raffan/impute_v2.3.2_x86_64_static/impute2";
$maps="/mnt/research2/anton/Eleanor_Raffan/maps";
$dog10k="/mnt/research2/anton/Eleanor_Raffan/dog10k/panel";
$rcmd="$impute -use_prephased_g -m $maps/chrXX.cf3.1_map.txt -h $dog10k/dog10k_plink_chrXX.phased.impute.haplotypes -l $dog10k/dog10k_plink_chrXX.phased.impute.legend -known_haps_g broad_plink_chrXX_gwas.phased.haps -int ZZ WW  -allow_large_regions -Ne 200 -o impute_chrXX/broad_plink_chrXX_gwas.phased.impute_final_chunkPP -phase";
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

for ($i=1;$i<=38;$i++){
	print "CHROMOSOME:$i CHRSIZE: $chr_size{$i}\n";
	$chunks=ceil($chr_size{$i}/$chunksize);
	print "Will divide into $chunks chunks\n";
	$size=ceil($chr_size{$i}/$chunks);

	if (-d "impute_chr$i"){
	} else {
		print "Making Dir: impute_chr$i\n";
		mkdir("impute_chr$i");
	}

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
