#!/usr/bin/perl
#
# extract the chromosome 8 sequences from the genome assemblies
#

open(IN, "Chrom8Scaffolds.txt") or die "failed to read the infile\n";

while(<IN>){
	chomp;
	@line = split(/\s+/,$_);
	open(FQ, "$line[2]") or die "failed to read $line[2]\n";
	open(OUT, "> chr8haplotypes/ch8_$line[0].fasta") or die "failed to write for $line[0]\n";
	$flag = 0;
	while(<FQ>){
		chomp;
		if(m/^>(\S+)/){
			$scaf = $1;
			if($scaf eq $line[1]){ ## chrom 8
				$flag = 1;
				print OUT ">$line[0]_ch8\n";
			} else{
				$flag = 0;
			}
		} elsif($flag==1){
			print OUT "$_\n";
		}
	}
	close(FQ);
	close(OUT);
}
