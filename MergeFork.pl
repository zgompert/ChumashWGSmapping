#!/usr/bin/perl
#
# merge alignments for each population sample with samtools version 1.16 
#


use Parallel::ForkManager;
my $max = 20;
my $pm = Parallel::ForkManager->new($max);

open(IDS,"pids.txt");
while(<IDS>){
	chomp;
	push(@IDs,$_);
}
close(IDS);

FILES:
foreach $id (@IDs){
	$pm->start and next FILES; ## fork
        system "samtools merge -c -p -o $id.bam $id"."_*.bam\n";
	system "samtools index -@ 2 $id.bam\n";
	$pm->finish;
}

$pm->wait_all_children;



