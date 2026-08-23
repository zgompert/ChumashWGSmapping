#!/usr/bin/perl
#
# syri comparative alignments
#

use Parallel::ForkManager;
my $max = 8;
my $pm = Parallel::ForkManager->new($max);

@rgenomes = ("ch8_24_0159_hap1.fasta","ch8_24_0163_hap1.fasta");


foreach $rg (@rgenomes){ 
	foreach $qg (@ARGV){
		$pm->start and next;
		$rg =~ m/ch8_([0-9a-z_]+)/;
		$rgid = $1;
		$qg =~ m/ch8_([0-9a-z_]+)/;
		$qgid = $1;
		$out = "syri_"."$rgid"."_"."$qgid";;
		## whole genome alignment with minimap2
		system "minimap2 -ax asm5 --eqx $rg $qg > $out.sam\n";
		## syri, -k keeps intermediate files, -F S is sam inpute
		system "syri -c $out.sam -r $rg -q $qg --prefix $out -F S --nosnp\n";
		## plot results from syri
		open(G,"> genomes$out") or die;
		print G "$rg\t$rgid\n";
		print G "$qg\t$qgid\n";
		close(G);
		system "plotsr --sr $out"."syri.out --genomes genomes$out -H 8 -W 5 -o $out"."_plot.pdf\n";

		$pm->finish;
	}
}
$pm->wait_all_children;

