#!/usr/bin/perl
#
# downsample.pl - Anton Enright (aje39@cam.ac.uk).
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

my $outfile=$ARGV[0];

$outfile=~ s/.vcf.gz/.downsampled.vcf/;
print "Processing $ARGV[0] to $outfile\n";

open(FILE,"../wisdom/markers_v5.csv");
while(<FILE>){
	chomp;
	my ($chr,$pos)=split("_",$_);
	$keep{$chr}{$pos}=1;

}


open(FILEOUT,">$outfile");

open(FILE,"gunzip -c $ARGV[0] |");
while(<FILE>){
	chomp;
	if (/^#/){
		print FILEOUT "$_\n";
	} else {
		@array=split("\t",$_);
		$keepsnp=0;
		if ($keep{$array[0]}{$array[1]}==1){
			print FILEOUT "$_\n";
		}
	}
}
close FILEOUT;

system("gzip $outfile");
