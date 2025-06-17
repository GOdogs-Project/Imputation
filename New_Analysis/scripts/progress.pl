#!/usr/bin/perl
#
# progress.pl - Anton Enright (aje39@cam.ac.uk).
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

print "|Chromosome\t|File\t|Size (MB)\t|Status\t|\n";
print "|-----------|-----|-----------|-------|\n";

open(PROC," du -sh dog10k_plink_chr*.refpanel.bed | sort -k1n|");
while(<PROC>){
	chomp;
	#230M	dog10k_plink_chr2.refpanel.bed
	@array=split("\t",$_);
	$size=$array[0];
	$file=$array[1];

	$chr=$file;
	$chr=~s/^.*chr(\d+).refpanel.*$/$1/g;
	$line=`tail -1 sbatch_shapeit_$chr.out`;
	chomp($line);

	if ($line=~ /Running time: (\d+) seconds/){
		$time=`grep \"Running time\" sbatch_shapeit_$chr.out | head -1`;
		chomp($time);
		$time=~ s/^Running time: (\d+) seconds/$1/g;
		$time=($1/60)/60;
		$time=sprintf("%2.2f",$time);
		$out="Finished: $time (hours)";
	} else {
		if ($line =~ /ERROR/){
			$line="";
		} else {
			#Burn-in iteration
			$other=$line;
			$other=~ s/\b//g;
			$other=~ s/[^A-Za-z0-9\[\]\/ ]//g;
			$other=~ s/\]([A-Za-z])/\|$1/g;
			@test=split(/\|/,$other);
			#print ">>$other<\n";
			$out="$test[-1]";
		}

	if ($line eq ''){
		$out="Not Started";
	}
	}

	$linebuffer{$chr}="|$chr\t|$file\t|$size\t|$out\t|\n";
}

foreach $thing(sort special(keys(%linebuffer))){
	print "$linebuffer{$thing}";
}

sub special{
	return($a<=>$b);
}
