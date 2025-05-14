#!/usr/bin/perl
#
# samtools/bcftools variant calling by LG 
#


use Parallel::ForkManager;
my $max = 10;
my $pm = Parallel::ForkManager->new($max);

my $genome ="/uufs/chpc.utah.edu/common/home/gompert-group4/data/timema/hic_genomes/t_chumash/hic/mod_timema_chumash_29Feb2020_N4ago.fasta";


foreach $chrom (@ARGV){
	$pm->start and next; ## fork
        $chrom =~ /chrom([0-9\.\-a-z]+)/ or die "failed here: $chrom\n";
	#$chrom =~ /chrom([0-9\.]+)/ or die "failed here: $chrom\n";
	$out = "o_tchum_chrom$1";
	system "bcftools mpileup -b bams -d 1000 -f $genome -R $chrom -a FORMAT/DP,FORMAT/AD -q 20 -Q 30 -I -Ou | bcftools call -v -c -p 0.01 -Ov -o $out"."vcf\n";


	$pm->finish;

}

$pm->wait_all_children;
