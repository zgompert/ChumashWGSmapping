#!/usr/bin/perl
#
# fit gemma BSLMM for T. chumash color 
#

use Parallel::ForkManager;
my $max = 4;
my $pm = Parallel::ForkManager->new($max);

$gfile1="hsub_mlpntest_tchum_red.geno";
$pfile1="ph_tchum_red.txt";
$gfile2="hsub_mlpntest_tchum_green.geno";
$pfile2="ph_tchum_green.txt";


foreach $ph (2){ 
	#foreach $ph (1..2){ 
	foreach $ch (20..59){
		#foreach $ch (0..59){
		sleep 2;
		$pm->start and next;
		## red genome
		$o = "o_poly_bslmm_red_ph$ph"."_ch$ch";
		if($ph == 1){
	    		print "gemma -g $gfile1 -p $pfile1 -bslmm 1 -n $ph -o $o -maf 0.01 -w 200000 -s 1000000 -k output/o_red_RG.cXX.txt\n";
	    		system "gemma -g $gfile1 -p $pfile1 -bslmm 1 -n $ph -o $o -maf 0.01 -w 200000 -s 1000000 -k output/o_red_RG.cXX.txt\n";
		} elsif($ph == 2){
	    		print "gemma -g $gfile1 -p $pfile1 -bslmm 1 -n $ph -o $o -maf 0.01 -w 200000 -s 1000000 -k output/o_red_GB.cXX.txt\n";
	    		system "gemma -g $gfile1 -p $pfile1 -bslmm 1 -n $ph -o $o -maf 0.01 -w 200000 -s 1000000 -k output/o_red_GB.cXX.txt\n";
		}
		## green genome
		$o = "o_poly_bslmm_green_ph$ph"."_ch$ch";
		if($ph == 1){
	    		print "gemma -g $gfile2 -p $pfile2 -bslmm 1 -n $ph -o $o -maf 0.01 -w 200000 -s 1000000 -k output/o_green_RG.cXX.txt\n";
	    		system "gemma -g $gfile2 -p $pfile2 -bslmm 1 -n $ph -o $o -maf 0.01 -w 200000 -s 1000000 -k output/o_green_RG.cXX.txt\n";
		} elsif($ph == 2){
	    		print "gemma -g $gfile2 -p $pfile2 -bslmm 1 -n $ph -o $o -maf 0.01 -w 200000 -s 1000000 -k output/o_green_GB.cXX.txt\n";
	    		system "gemma -g $gfile2 -p $pfile2 -bslmm 1 -n $ph -o $o -maf 0.01 -w 200000 -s 1000000 -k output/o_green_GB.cXX.txt\n";
		}
		$pm->finish;
	}
}
$pm->wait_all_children;

