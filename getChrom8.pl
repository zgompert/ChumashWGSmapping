#!/usr/bin/perl
#
# extract chromosome 8 from T. chumash genomes
# for all ten haplotypes, this is scaffold 4 
#
## give it: 24_01*/Hap*Chr.fasta or equivalent masked
## I am trying both ways

foreach $fa (@ARGV){
	print "Working on $fa\n";
	$fa =~ m!([0-9_]+)/Hap(\d)Chr\.([a-z\.]+)! or die;
	$out = "ch8_$1_hap$2.$3\n";
	open(IN, $fa) or die;
	open(OUT, "> $out") or die "failed to write\n";
	$flag = 0;
	while(<IN>){
		chomp;
		if(s/^>hap\d_Chr04/>Chr8/){ ## scaffold 4 = chromosome 8
			$flag = 1;
		} elsif	(m/^>/){ ## some other scaffold
			$flag = 0;
		}
		if($flag == 1){
			print OUT "$_\n";
		}
	}
	close(IN);
	close(OUT);
}
