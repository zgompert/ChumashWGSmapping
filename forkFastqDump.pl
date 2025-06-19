#!/usr/bin/perl
#
# sra to fastq 
#


use Parallel::ForkManager;
my $max = 24;
my $pm = Parallel::ForkManager->new($max);

DIR:
foreach $sra (@ARGV){
	$pm->start and next DIR; ## fork
        system "fasterq-dump ./$sra\n";
	$pm->finish;
}

$pm->wait_all_children;



