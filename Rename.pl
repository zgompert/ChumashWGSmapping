#!/usr/bin/perl
#
# simplify names
#

foreach $fq (@ARGV){
	$fq =~ m/19_(\d+)_R([12])\.fastq$/ or die "failed here: $fq\n";
	$a = $1;
	$b = $2;
	$o = "tchum_19_$a"."_$b.fastq";
	system "mv $fq $o\n";
}
