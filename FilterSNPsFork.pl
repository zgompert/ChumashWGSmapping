#!/usr/bin/perl

use Parallel::ForkManager;
my $max = 4;
my $pm = Parallel::ForkManager->new($max);

FILES:
foreach $vcf (@ARGV){
        $pm->start and next FILES; ## fork
	$vcf =~ m/(\d)/;
	$ch = $1;
	system "bgzip -d $vcf\n";
	$snpf = "comb_filt_ch$ch.txt";
	open(IN,$snpf) or die "failed to read snpfile: $snpf\n";
	while(<IN>){
		chomp;
		m/(^hap\S+\s+\S+)/;
		$drop{$1} = 1;
	}
	close(IN);
	$vcf =~ s/\.gz//;
	open(IN,$vcf) or die "failed to read the vcf: $vcf\n";
	open(OUT, "> clean_$vcf") or die "failed ot write to clean_$vcf\n";
	$cnt =0;
	$kep =0;
	while(<IN>){
		if(m/^#/){
			print OUT $_;
		} else{
			m/(^hap\S+\s+\S+)/ or print "warning, what is this $_\n";
			$s = $1;
			$cnt++;
			unless(defined $drop{$s}){
				print OUT $_;
				$kep++;
			}
		}
	}
	close(IN);
	close(OUT);
	print "for $vcf, kept $kep out of $cnt SNPs\n";
	print "bgzip -c clean_$vcf\n";
        $pm->finish;
}

$pm->wait_all_children;

