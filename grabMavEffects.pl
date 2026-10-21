#!/usr/bin/perl

$ordfile = shift(@ARGV);
open(OSNP, $ordfile) or die;
while(<OSNP>){
        chomp;
        push(@osnps,$_);
}

$snpfile = shift(@ARGV);

open(SNP, $snpfile) or die "could not open the SNP file\n";
$pos = 0;
while(<SNP>){
        chomp;
        @line = split(/\s+/,$_);
        if($line[1] == $osnps[$pos]){
                push(@snps,$_);
                $pos++;
        }
}
close(SNP);
print "Ended on SNP number $pos\n";

foreach $in (@ARGV){
	$pos = 0;
	foreach $snp (@snps){ ## define and 0 all SNP effects
		$eff{$pos} = 0;
		$pos++
	}
	## combine information over reps and sort
	open(OUT, "> mav_$in") or die "failed to write mav $in\n";
	$nreps=60;
	foreach $rep (0..59){
		$in =~ s/ch\d+/ch$rep/ or die "failed sub for $in to $rep\n";
		print "working on $in\n";
		open(IN, $in) or die "failed to open the infile $in\n";
		<IN>; ## burn header
		$pos = 0;
		while(<IN>){
			chomp;
			@line = split(/\s+/,$_);
			$eff{$pos} += ($line[4] + $line[5] * $line[6]) * (1/$nreps);
			$pos++;
		}
		close(IN);
	}
	$pos = 0;
	foreach $snp (@snps){
		print OUT "$snp $eff{$pos}\n";
		$pos++;
	}
	close(OUT);
}
