#!/usr/bin/perl
#
# liftover.pl - Anton Enright (aje39@cam.ac.uk).
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

$cmd="../liftover/bcftools-1.22/bcftools +liftover --no-version -Ou XXXX -- -c ../liftover/canFam3ToCanFam4.over.chain.gz -s ../liftover/canFam3.fa -f ../liftover/canFam4.fa | ../liftover/bcftools-1.22/bcftools sort -Oz -o YYYY";

foreach $file (@ARGV){
	print "FILE: $file\n";
	$stem=$file;
	$stem=~ s/.vcf.gz//g;
	print "STEM: $stem\n";
	$outfile="/tmp/$stem.cf4.vcf.gz";
	open(FILEOUT,">/tmp/$stem.vcf");
	open(PROC,"gunzip -c $file |");
	while(<PROC>){
		chomp;
		$line=$_;
		if (/\#\#contig\=<ID=(\d+)/){
			$chr=$1;
			$line=~ s/ID=$chr/ID=chr$chr/g;
		}

		if (/^(\d+)\s+/){
			$line=~ s/^(\d+)/chr$1/g;
		}

		print FILEOUT "$line\n";
	}
	close(FILEOUT);


	$rcmd=$cmd;
	$rcmd=~ s/XXXX/\/tmp\/$stem.vcf/g;
	$rcmd=~ s/YYYY/$outfile/g;

	print "Running: $rcmd\n";
	system("$rcmd");

	open(FILEOUT,"> /tmp/$stem.cf4.vcf");
	open(PROC,"gunzip -c $outfile|");
	while(<PROC>){
		chomp;
		$line=$_;
		if (/\#\#contig\=<ID=chr(\d+)/){
			$line=~ s/ID=chr/ID=/g;
		}
		if (/^chr(\d+)\s+/){
			$line=~ s/^chr(\d+)/$1/g;
		}
		print FILEOUT "$line\n";

	}
	close FILEOUT;
	print "Zipping: /tmp/$stem.cf4.vcf\n";
	print "Moving: /tmp/$stem.cf4.vcf.gz to .\n";
	system("gzip --force /tmp/$stem.cf4.vcf; mv /tmp/$stem.cf4.vcf.gz .");
}
