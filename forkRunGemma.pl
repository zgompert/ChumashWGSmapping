#!/usr/bin/perl
#
# fit gemma BSLMM for T. chumash color 
#

use Parallel::ForkManager;
my $max = 40;
my $pm = Parallel::ForkManager->new($max);

$g = "hsub_mlpntest_EXP_tchum.geno";
$p = "ph_tchum2019.txt";


foreach $ph (1..2){ 
	foreach $ch (0..29){
		sleep 2;
		$pm->start and next;
		$o = "o_poly_bslmm_ph$ph"."_ch$ch";
		if($ph == 1){
	    		system "gemma -g $g -p $p -bslmm 1 -n $ph -o $o -maf 0 -w 200000 -s 1000000 -k output/o_poly_RG.cXX.txt\n";
		} elsif($ph == 2){
	    		system "gemma -g $g -p $p -bslmm 1 -n $ph -o $o -maf 0 -w 200000 -s 1000000 -k output/o_poly_GB.cXX.txt\n";
		}
		$pm->finish;
	}
}
$pm->wait_all_children;

