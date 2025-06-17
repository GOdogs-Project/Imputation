#!/usr/bin/perl
#
# shape_it.pl - Anton Enright (aje39@cam.ac.uk).
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
$run=0;
$batch="sbatch --partition=PlanEx2 --mem=60G --cpus-per-task 20 --ntasks 1 -o XX -e YY";
foreach $arg(@ARGV){
	if ($arg eq '--runnit'){
		$run=1;
	}
}

$cmd="stdbuf -oL -eL srun -u ../../shapeit.v2.904.3.10.0-693.11.6.el7.x86_64/bin/shapeit -M ../../maps/chrXX.cf3.1_map.txt -B dog10k_plink_chrXX.refpanel -O dog10k_plink_chrXX.phased -T 20 --window 2 --effective-size 200 --force\nstdbuf -oL -eL srun -u ../../shapeit.v2.904.3.10.0-693.11.6.el7.x86_64/bin/shapeit -convert --input-haps dog10k_plink_chrXX.phased --output-ref dog10k_plink_chrXX.phased.impute";
for ($i=1;$i<=38;$i++){

	$rcmd=$cmd;
	$rcmd=~ s/XX/$i/g;
	print "$rcmd\n";
	open(FILE,">/tmp/$$\_sbatch$i.sh");
        print FILE "#!/bin/sh\n$rcmd\n";
	$sbatch=$batch;
	$sbatch=~ s/XX/sbatch_shapeit_$i.out/g;
	$sbatch=~ s/YY/sbatch_shapeit_$i.err/g;
	print "$sbatch /tmp/$$\_sbatch$i.sh\n";
	if ($run){
		system("$sbatch /tmp/$$\_sbatch$i.sh");
	}
}

