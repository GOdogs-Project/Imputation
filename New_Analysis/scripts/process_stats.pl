#!/usr/bin/perl
#
# process_stats.pl - Anton Enright (aje39@cam.ac.uk).
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

my @r2=();
my @con=();
my @info=();

print "chrom\tchunk\tr2_med\tr2_avg\tr2_std\tcon_med\tcon_avg\tcon_std\tinf_med\tinf_avg\tinf_std\n";
foreach $file(@ARGV){
	open(FILE,$file);
	my @r2_chunk=();
	my @con_chunk=();
	my @info_chunk=();

	$chunk_no=$file;
	$chunk_no=~ s/.*chunk(\d+).*/$1/g;
	$chr=$file;
	$chr=~ s/.*chr(\d+).*/$1/g;
	while(<FILE>){
		chomp;
		@array=split(" ",$_);
		if ($array[0] ne 'snp_id'){
			if (($array[1] ne ".") && ($array[11] != -1)){
				push(@r2_chunk,$array[11]);
				push(@con_chunk,$array[10]);
				push(@info_chunk,$array[9]);
			}
		}

	}
	push(@r2,@r2_chunk);
	push(@con,@con_chunk);
	push(@info,@info_chunk);
	print "$chr\t$chunk_no\t";
	print median(@r2_chunk) . "\t" . sprintf("%2.4f",get_avg(\@r2_chunk)) ."\t" . get_stddev(\@r2_chunk) . "\t";
	print median(@con_chunk) . "\t" . sprintf("%2.4f",get_avg(\@con_chunk)) ."\t" .get_stddev(\@con_chunk) . "\t";
	print median(@info_chunk) . "\t".  sprintf("%2.4f",get_avg(\@info_chunk)) ."\t" .get_stddev(\@info_chunk) . "\n";
}

print "$chr\tTOTAL\t";
print median(@r2) . "\t" . sprintf("%2.4f",get_avg(\@r2)) ."\t" . get_stddev(\@r2) . "\t";
print median(@con) . "\t" . sprintf("%2.4f",get_avg(\@con)) ."\t" .get_stddev(\@con) . "\t";
print median(@info) . "\t".  sprintf("%2.4f",get_avg(\@info)) ."\t" .get_stddev(\@info) . "\n";


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
