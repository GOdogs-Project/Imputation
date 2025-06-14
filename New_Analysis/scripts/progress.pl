#!/usr/bin/perl

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
