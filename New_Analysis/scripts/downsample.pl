#!/usr/bin/perl

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
