#!/usr/bin/perl
#
# process_stats_individual.pl - Anton Enright (aje39@cam.ac.uk).
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

$iterator=0;
#cat broad_plink_chr31_gwas.phased.impute.samples | cut -f 2 -d ' ' |grep -v 'ID_2'

for ($i=1;$i<=38;$i++){
	open(PROC,"cat broad_plink_chr$i\_gwas.phased.impute.samples | cut -f 2 -d \' \' | grep -v \'ID_2\'|");
	$iterator=0;
	while(<PROC>){
		chomp;
		$sample_name{$i}{$iterator}=$_;
		$iterator++;
		$chromosomes{$_}++;
	}

}
#open(FILE,"samples_imputed.txt");
#while(<FILE>){
#chomp;
#$sample_name{$iterator}=$_;
#$iterator++;
#}

open(FILE,"sample_info.tsv");
while(<FILE>){
chomp;

$line=$_;
$line=~ s/[\"\']//g;
@array=split("\t",$line);
$field1{$array[0]}=$array[1];
$field2{$array[0]}=$array[2];
$field3{$array[0]}=$array[3];

}

foreach $dir (@ARGV){
	undef(%concord);
	undef(%r2);
foreach $file(`ls $dir/*_info_by_sample`){
	chomp($file);
	open(FILE,$file);
	$chunk_no=$file;
	$chunk_no=~ s/.*chunk(\d+).*/$1/g;
	$chr=$file;
	$chr=~ s/.*chr(\d+).*/$1/g;

	$individual=0;
	while(<FILE>){
		chomp;
		@array=split(" ",$_);
		if($array[0] ne 'concord_type0'){
			#print ">$_\n";
			#print "[$individual] [$sample_name{$chr}{$individual}] [$field2{$sample_name{$chr}{$individual}}]@array\n";
			push(@{$concord{$individual}},$array[0]);
			push(@{$r2{$individual}},$array[1]);
			$individual++;
		}
	}
}

for ($i=0;$i<$individual;$i++){
	print "$chr\t$i\t$sample_name{$chr}{$i}\t$field1{$sample_name{$chr}{$i}}\t$field2{$sample_name{$chr}{$i}}\t$field3{$sample_name{$chr}{$i}}\t";
	print median(@{$r2{$i}}) . "\t";
	print median(@{$concord{$i}}) . "";
	print "\n";
}

}


sub get_stddev {
  return sprintf("%2.4f",sqrt(get_disp(@_)));
}

sub get_disp {
  my ($array_ref) = @_;
  my $mean = get_avg($array_ref);
  my $count = @$array_ref;
  my $sum = 0;
  
  for my $num (@$array_ref) {
      $sum += (($num - $mean) ** 2);
  }
  return $sum / $count;
}

sub get_avg {
  my ($array_ref) = @_;
  my $count = @$array_ref;
  my $sum = 0;
  for my $num (@$array_ref) {
      $sum += $num;
  }
  return $sum / $count;
}

sub median
{
    my @vals = sort {$a <=> $b} @_;
    my $len = @vals;
    if($len%2) #odd?
    {
        return $vals[int($len/2)];
    }
    else #even
    {
        return (sprintf("%2.4f",$vals[int($len/2)-1] + $vals[int($len/2)])/2);
    }
}
