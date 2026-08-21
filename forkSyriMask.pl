#!/usr/bin/perl
#
# syri comparative alignments, using masked genomes
#

use Parallel::ForkManager;
my $max = 8;
my $pm = Parallel::ForkManager->new($max);

@rgenomes = ("ch8_24_0159_hap1.fasta.masked","ch8_24_0163_hap1.fasta.masked");


foreach $rg (@rgenomes){ 
	foreach $qg (@ARGV){
		$pm->start and next;
		$rg =~ m/ch8_([0-9a-z_]+)/;
		$rgid = $1;
		$qg =~ m/ch8_([0-9a-z_]+)/;
		$qgid = $1;
		$out = "syri_"."$rgid"."_"."$qgid"."_mask";
		## whole genome alignment with minimap2
		system "minimap2 -ax asm5 --eqx $rg $qg > $out.sam\n";
		## syri, -k keeps intermediate files, -F S is sam inpute
		system "syri -c $out.sam -r $rg -q $qg -k -F S --nosnp\n";
		## plot results from syri
		system "plotsr $out"."_syri.out $rg $qg -H 8 -W 5\n";

		$pm->finish;
	}
}
$pm->wait_all_children;

