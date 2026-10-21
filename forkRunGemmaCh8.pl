#!/usr/bin/perl
#
# fit gemma BSLMM for T. chumash color 
# ch8 only

use Parallel::ForkManager;
my $max = 8;
my $pm = Parallel::ForkManager->new($max);

$gfile1="ch8_mlpntest_tchum_red.geno";
$pfile1="ph_tchum_red.txt";
$gfile2="ch8_mlpntest_tchum_green.geno";
$pfile2="ph_tchum_green.txt";

## kinship matrix first
system "gemma -g  $gfile1 -p $pfile1 -gk 1 -o ch8_red -maf 0.01\n";
#system "gemma -g  $gfile2 -p $pfile2 -gk 1 -o ch8_green -maf 0.01\n";

foreach $ph (1..2){ 
	foreach $ch (0..19){
		sleep 2;
		$pm->start and next;
		## red genome
		$o = "o_ch8_poly_bslmm_red_ph$ph"."_ch$ch";
		if($ph == 1){
	    		system "gemma -g $gfile1 -p $pfile1 -bslmm 1 -n $ph -o $o -maf 0.01 -w 200000 -s 1000000 -k output/ch8_red.cXX.txt\n";
		} elsif($ph == 2){
	    		system "gemma -g $gfile1 -p $pfile1 -bslmm 1 -n $ph -o $o -maf 0.01 -w 200000 -s 1000000 -k output/ch8_red.cXX.txt\n";
		}
		## green genome
		#		$o = "o_ch8_poly_bslmm_green_ph$ph"."_ch$ch";
		#		if($ph == 1){
		#	    		system "gemma -g $gfile2 -p $pfile2 -bslmm 1 -n $ph -o $o -maf 0.01 -w 200000 -s 1000000 -k output/ch8_green.cXX.txt\n";
		#		} elsif($ph == 2){
		#	    		system "gemma -g $gfile2 -p $pfile2 -bslmm 1 -n $ph -o $o -maf 0.01 -w 200000 -s 1000000 -k output/ch8_green.cXX.txt\n";
		#		}
		$pm->finish;
	}
}
$pm->wait_all_children;

