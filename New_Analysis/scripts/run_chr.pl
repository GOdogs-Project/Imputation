#!/usr/bin/perl
#
# run_chr.pl - Anton Enright (aje39@cam.ac.uk).
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

