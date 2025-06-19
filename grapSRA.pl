#!/usr/bin/perl
#

# ml sra-toolkit

open(IN,"SRAids2.txt") or die;
#open(IN,"SRAids.txt") or die;
while(<IN>){
	chomp;
	system "prefetch $_\n";
}
