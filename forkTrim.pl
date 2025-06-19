#!/usr/bin/perl
#
# filter RNA sequences with trimmomatic 
#


use Parallel::ForkManager;
my $max = 24;
my $pm = Parallel::ForkManager->new($max);

FILE:
foreach $fq1 (@ARGV){
	$pm->start and next FILE; ## fork
	$fq2 = $fq1;
	$fq2 =~ s/_1\.fastq/_2.fastq/ or die "failed here\n";
	print "trimmomatic PE -threads 3 -phred33 $fq1 $fq2 clean_$fq1 unpaired_$fq1 clean_$fq2 unpaired_$fq2 ILLUMINACLIP:Illumina_TruSeq.fa:2:30:8 LEADING:20 TRAILING:20 SLIDINGWINDOW:4:20 MINLEN:35\n";
	system "trimmomatic PE -threads 3 -phred33 $fq1 $fq2 clean_$fq1 unpaired_$fq1 clean_$fq2 unpaired_$fq2 ILLUMINACLIP:Illumina_TruSeq.fa:2:30:8 LEADING:20 TRAILING:20 SLIDINGWINDOW:4:20 MINLEN:35\n";
	#print "trimmomatic PE -threads 3 -phred33 $fq1 $fq2 clean_$fq1 unpaired_$fq1 clean_$fq2 unpaired_$fq2 ILLUMINACLIP:Illumina_TruSeq.fa:2:30:8:1:true LEADING:20 TRAILING:20 SLIDINGWINDOW:4:20 MINLEN:35\n";
	
	$pm->finish;
}
## removes leading or trailing low quality bases below qual 20
### drops reads < 35 bps
### cut when 4 bp window has quality less than 20 on average
### adapter uses the specified file, allows for 2 missmatches to try full match, palindrome clip of 30, and simple clip of 8
$pm->wait_all_children;



