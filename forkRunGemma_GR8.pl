#!/usr/bin/perl
#
# fit gemma BSLMM for T. chumash color 
#

use Parallel::ForkManager;
my $max = 4;
my $pm = Parallel::ForkManager->new($max);

$gfile2="hsub_mlpntest_GR8_tchum_gr8_green.geno";
$pfile2="ph_tchumGR8.txt";


foreach $ph (1..2){ 
	foreach $ch (0..59){
		sleep 2;
		$pm->start and next;
		## green genome
		$o = "o_poly_bslmm_GR8_ph$ph"."_ch$ch";
		if($ph == 1){
	    		print "gemma -g $gfile2 -p $pfile2 -bslmm 1 -n $ph -o $o -maf 0.01 -w 200000 -s 1000000 -k output/o_GR8_RG.cXX.txt\n";
	    		system "gemma -g $gfile2 -p $pfile2 -bslmm 1 -n $ph -o $o -maf 0.01 -w 200000 -s 1000000 -k output/o_GR8_RG.cXX.txt\n";
		} elsif($ph == 2){
	    		print "gemma -g $gfile2 -p $pfile2 -bslmm 1 -n $ph -o $o -maf 0.01 -w 200000 -s 1000000 -k output/o_GR8_GB.cXX.txt\n";
	    		system "gemma -g $gfile2 -p $pfile2 -bslmm 1 -n $ph -o $o -maf 0.01 -w 200000 -s 1000000 -k output/o_GR8_GB.cXX.txt\n";
		}
		$pm->finish;
	}
}
$pm->wait_all_children;

