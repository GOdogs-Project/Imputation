#!/usr/bin/perl
#
# process_persnp_stats.pl - Anton Enright (aje39@cam.ac.uk).
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

open(PROC,"cat impute_chr*/broad_plink_chr*_gwas.phased.impute_final_chunk*_info | grep -v \'\\-\\-\\-\' | grep -v \'\\-1\' | cut -f 1,3,6,10,11,12 -d \' \'|");
open(FILEOUT,">all_stats_by_chr.txt");
while(<PROC>){
	chomp;
	#snp_id position exp_freq_a1 info_type0 concord_type0 r2_type0
	#10 452549 0.604 0.892 0.951 0.914
	#10 523567 0.710 0.955 0.973 0.942
	#10 670696 0.701 0.958 0.978 0.945

	if (!(/^snp_id/)){
	$line=$_;
	$line=~ s/ /\t/g;
	print FILEOUT "$line\n";
	#$snp_hash_r2{$array[0]}{$array[1]}=$array[5];
	#$snp_hash_concord{$array[0]}{$array[1]}=$array[4];
	#$snp_hash_info{$array[0]}{$array[1]}=$array[3];
	}
}

close(FILEOUT);

sub numeric{
	return($a <=> $b);
}
