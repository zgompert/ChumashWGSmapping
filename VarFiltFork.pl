#!/usr/bin/perl
#
# filter vcf files 
#


use Parallel::ForkManager;
my $max = 48;
my $pm = Parallel::ForkManager->new($max);


foreach $vcf (@ARGV){
	$pm->start and next; ## fork
	print "perl vcfFilter.pl $vcf\n";
	system "perl vcfFilter.pl $vcf\n";
	$pm->finish;

}

$pm->wait_all_children;



