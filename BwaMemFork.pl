#!/usr/bin/perl
#
# alignment with bwa mem 
#


use Parallel::ForkManager;
my $max = 12;
my $pm = Parallel::ForkManager->new($max);
my $genome = "/uufs/chpc.utah.edu/common/home/gompert-group4/data/timema/hic_genomes/t_chumash/hic/mod_timema_chumash_29Feb2020_N4ago.fasta";

$files = shift(@ARGV);

open(IN, $files);
FILES:
while(<IN>){
	$pm->start and next FILES; ## fork
	chomp;
	$fq1 = $_;
	$fq2 = $fq1;
	$fq2 =~ s/R1\.fastq\.gz/R2.fastq.gz/ or die "failed substitution for $fq1\n";
        $fq1 =~ m/(^[a-zA-Z0-9\-]+)_/ or die "failed to match id $fq1\n";
	$ind = $1;
	$fq1 =~ m/([A-Za-z0-9]+_L\d+)_R1.fastq.gz$/ or die "failed to match library for $fq1\n";
	$lib = $1;
        print "/uufs/chpc.utah.edu/common/home/u6000989/source/bwa-mem2-2.0pre2_x64-linux/bwa-mem2 mem -t 1 -k 19 -r 1.5 -R \'\@RG\\tID:wgs-"."$ind\\tLB:wgs-"."$ind"."_$lib\\tSM:wgs-"."$ind"."\' $genome $fq1 $fq2 | samtools sort -@ 2 -O BAM -o $ind"."_$lib.bam - && samtools index -@ 2 $ind"."_$lib.bam\n";
        system "/uufs/chpc.utah.edu/common/home/u6000989/source/bwa-mem2-2.0pre2_x64-linux/bwa-mem2 mem -t 1 -k 19 -r 1.5 -R \'\@RG\\tID:wgs-"."$ind\\tLB:wgs-"."$ind"."_$lib\\tSM:wgs-"."$ind"."\' $genome $fq1 $fq2 | samtools sort -@ 2 -O BAM -o $ind"."_$lib.bam - && samtools index -@ 2 $ind"."_$lib.bam\n";

	$pm->finish;
}
close(IN);

$pm->wait_all_children;


