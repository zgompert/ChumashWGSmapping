# ChumashWGSmapping
Genomic analyses of color and survival in a selection experiment for *Timema chumash* based on whole genome sequence data

# Data

The DNA sequence data for this project are in two places. The first two batches are in `/uufs/chpc.utah.edu/common/home/gompert-group4/data/t_chumash_wgs/`: batch 1 = 192 samples, batch 2 = 288 samples (some of the latter are not part of the mapping and selection experiment). See [Nosil et al. 2020](https://www.proquest.com/docview/2473271380?pq-origsite=gscholar&fromopenview=true&sourcetype=Scholarly%20Journals) ([Nosil2020.pdf](https://github.com/user-attachments/files/17368867/Nosil2020.pdf)) for our analyses of GBS data from this same experiment. The data from the 3rd batch are in `/uufs/chpc.utah.edu/common/home/gompert-group5/data/TchumJan2025/`. Batch 3 includes 541 additioanl *T. chumash*

**batch 1:** These samples were sequenced by GeT (Genotoul) and delivered July 2024. The IDs are at the start of the file names (Tchum-XXX, where XXX denotes the digits after 19_ in the actual file names). There are two files each of forward and reverse reads for each stick insect. 

**batch 2:** These samples were sequenced by GeT (Genotoul) and delivered October 2024. The IDs are again at the start of the file name but here are either 19_XXX (the actual ID we use) or GR806-XX for the non-experimental stick insects.

**batch 3:** These samples were sequenced by GeT (Genotoul) and delivered January 2025. The IDs are again at the start of the file name but here are either 15_XXX. These are all non-experimental *T. chumash* from natural populations. The population information is in `/uufs/chpc.utah.edu/common/home/gompert-group5/data/TchumJan2025/sample_sp_loc_host.dsv'. The GR806-XX samples from batch 2 go with these individuals (i.e., they add to this data set).

The experimental data (color and survival) are in the 2019_Tchumash* files in `/uufs/chpc.utah.edu/common/home/gompert-group4/data/t_chumash_wgs/`.

I am using our 2020 *T. chumash* genome assembly. This genome is described in [Nosil et al. 2020](https://www.proquest.com/docview/2473271380?pq-origsite=gscholar&fromopenview=true&sourcetype=Scholarly%20Journals). See [timema_chumash_29Feb2020_N4ago.report.pdf](https://github.com/user-attachments/files/17368853/timema_chumash_29Feb2020_N4ago.report.pdf) for a report from Dovetail. Past comparative alignments suggest general chromosome-level synteny with *T. crstinae* except that *T. cristinae* chromsomes 1, 3, 4 and 9 appear to be fused in *T. chumash*. Patrik suggested this is consistent with cytology (add source). Also note (from the report) that the genome size is notably larger than our other *Timema* genomes (2.4 Gbps) and the final assembly has a non-trivial number of duplicated BUSCOs (29). This could reflect some difficulities in the assembly (especially if the two genome copies differed enough) that might be cleaned up with phased genomes later. However, chromosome 8 appears to be largely syntenic (though not colinear) between *T. cristinae* and *T. chumash* and is the main focus of our interest ([see the TimemaFusion repository for past work on this](https://github.com/zgompert/TimemaFusion)).

# DNA sequence alignment

I aligned the whole genome sequence data to the *T. chumash* genome using `bwa-mem2` (version 2.0pre2).

I first indexed the genome.

```bash
/uufs/chpc.utah.edu/common/home/u6000989/source/bwa-mem2-2.0pre2_x64-linux/bwa-mem2 index mod_timema_chumash_29Feb2020_N4ago.fasta
```
I then ran the alignments. I copied all of the fastq files to scratch space for this.

```bash
#!/bin/bash
#SBATCH --time=240:00:00
#SBATCH --nodes=1
#SBATCH --ntasks=20
#SBATCH --account=gompert-np
#SBATCH --partition=gompert-np
#SBATCH --job-name=bwa-mem2
#SBATCH --mail-type=FAIL
#SBATCH --mail-user=zach.gompert@usu.edu

module load samtools
##Version: 1.16 (using htslib 1.16)

cd /scratch/general/nfs1/u6000989/t_chumash_wgs

perl BwaMemFork.pl ./*R1.fastq.gz
```

This runs the following, which includes piping the results on to `samtools` (version 1.16) to compress, sort and index the alignments. These steps were done in a few batches.

```perl
#!/usr/bin/perl
#
# alignment with bwa mem 
#


use Parallel::ForkManager;
my $max = 40;
my $pm = Parallel::ForkManager->new($max);
my $genome = "/uufs/chpc.utah.edu/common/home/gompert-group4/data/timema/hic_genomes/t_chumash/hic/mod_timema_chumash_29Feb2020_N4ago.fasta";

FILES:
foreach $fq1 (@ARGV){
	$pm->start and next FILES; ## fork
	$fq2 = $fq1;
	$fq2 =~ s/R1\.fastq\.gz/R2.fastq.gz/ or die "failed substitution for $fq1\n";
        $fq1 =~ m/(^[a-zA-Z0-9\-]+)_/ or die "failed to match id $fq1\n";
	$ind = $1;
	$fq1 =~ m/([A-Za-z0-9]+_L\d+)_R1.fastq.gz$/ or die "failed to match library for $fq1\n";
	$lib = $1;
        system "/uufs/chpc.utah.edu/common/home/u6000989/source/bwa-mem2-2.0pre2_x64-linux/bwa-mem2 mem -t 1 -k 19 -r 1.5 -R \'\@RG\\tID:wgs-"."$ind\\tLB:wgs-"."$ind"."_$lib\\tSM:wgs-"."$ind"."\' $genome $fq1 $fq2 | samtools sort -@ 2 -O BAM -o $ind"."_$lib.bam - && samtools index -@ 2 $ind"."_$lib.bam\n";

	$pm->finish;
}
```
Next, I merge the bam files for each individual (when run in multiple lanes) to generate one bam per individual. I used `samtools` for this.

```bash
#!/bin/bash
#SBATCH --time=240:00:00
#SBATCH --nodes=1
#SBATCH --ntasks=20
#SBATCH --account=gompert-np
#SBATCH --partition=gompert-np
#SBATCH --job-name=merge
#SBATCH --mail-type=FAIL
#SBATCH --mail-user=zach.gompert@usu.edu

module load samtools
##Version: 1.16 (using htslib 1.16)
cd /scratch/general/nfs1/u6000989/t_chumash_wgs

perl MergeFork.pl 
```
```perl

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
```
I then used `samtools` to mark and remove (`-r` option) PCR duplicates.

```bash
#!/bin/sh
#SBATCH --time=240:00:00
#SBATCH --nodes=1
#SBATCH --ntasks=12
#SBATCH --account=gompert-np
#SBATCH --partition=gompert-np
#SBATCH --job-name=dedup
#SBATCH --mail-type=FAIL
#SBATCH --mail-user=zach.gompert@usu.edu

module load samtools
##Version: 1.16 (using htslib 1.16)

cd /scratch/general/nfs1/u6000989/t_chumash_wgs


perl RemoveDupsFork.pl *bam
```
```perl
#!/usr/bin/perl
#
# PCR duplicate removal with samtools
#


use Parallel::ForkManager;
my $max = 24;
#my $max = 40;
my $pm = Parallel::ForkManager->new($max);

FILES:
foreach $bam (@ARGV){
	$pm->start and next FILES; ## fork
	$bam =~ m/^([A-Za-z0-9_\-]+)/ or die "failed to match $bam\n";
	$base = $1;
	system "samtools collate -o co_$base.bam $bam\n";
	system "samtools fixmate -m co_$base.bam fix_$base.bam\n";
	system "samtools sort -o sort_$base.bam fix_$base.bam\n";
	## using default definition of dups
	## measure positions based on template start/end (default). = -m t
	system "samtools markdup -T /scratch/general/nfs1/u6000989/t_chumash_wgs/Dedup/ -r sort_$base.bam dedup_$base.bam\n";
	$pm->finish;
}

$pm->wait_all_children;
```

This left me with a total of 1020 bam files (one per individual), which can be found in `/uufs/chpc.utah.edu/common/home/gompert-group5/data/t_chumash_wgs/` (batches 1 and 2) and `/uufs/chpc.utah.edu/common/home/gompert-group5/data/TchumJan2025/` (batch). From this point on, both batches were processed together.

# Variant calling and filtering

I used `bcftools` (version 1.16) to identifty SNPs. This was done in parallel across chromsomes (with chromsome 1 split into four subsets for more efficient processing). The chrom*list files give the chromosome IDs and regions for chromosome 1. The bams file lists all 1020 (deduplicated) bams.

```bash
#!/bin/sh
#SBATCH --time=240:00:00
#SBATCH --nodes=1
#SBATCH --ntasks=24
#SBATCH --account=gompert-np
#SBATCH --partition=gompert-np
#SBATCH --job-name=bcf_call
#SBATCH --mail-type=FAIL
#SBATCH --mail-user=zach.gompert@usu.edu


module load samtools
## version 1.16
module load bcftools
## version 1.16


cd /scratch/general/nfs1/u6000989/t_chumash_wgs

perl BcfForkLg.pl chrom*list 
```
```perl
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
```

I then filtered the SNP set for coverage, missing data, and various tests of bias using my own perl script.

```perl
#!/usr/bin/perl

use warnings;
use strict;

# this program filters a vcf file based on overall sequence coverage, number of non-reference reads, number of alleles, and reverse orientation reads

# usage vcfFilter.pl infile.vcf
#
# change the marked variables below to adjust settings
#

#### stringency variables, edits as desired
## 1020 inds, 2x
my $minCoverage = 2040; # minimum number of sequences; DP
my $minAltRds = 10; # minimum number of sequences with the alternative allele; AC
my $notFixed = 1.0; # removes loci fixed for alt; AF
my $bqrs = 3; # Z-score base quality rank sum test; BaseQRankSum
my $mqrs = 3; # Z-score mapping quality rank sum test; MQRankSum
my $rprs = 3; # Z-score read position rank sum test; ReadPosRankSum
my $mq = 30; # minimum mapping quality; MQ
my $miss = 102; # maximum number of individuals with no data = 10%
##### this set is for GBS
my $d;

my @line;

my $in = shift(@ARGV);
open (IN, $in) or die "Could not read the infile = $in\n";
$in =~ m/^([a-zA-Z_0-9\-]+)\.listvcf$/ or die "Failed to match the variant file\n";
open (OUT, "> filtered2x_$1.vcf") or die "Could not write the outfile\n";

my $flag = 0;
my $cnt = 0;

while (<IN>){
	chomp;
	$flag = 1;
	if (m/^\#/){ ## header row, always write
		$flag = 1;
	}
	elsif (m/^Sc/){ ## this is a sequence line, you migh need to edit this reg. expr.
		$flag = 1;
		$d = () = (m/\d\/\d:0,0,0:0/g); ## for bcftools call
		if ($d >= $miss){
			$flag = 0;
			##print "fail missing : ";
		}
		if (m/[ACTGN]\,[ACTGN]/){ ## two alternative alleles identified
			$flag = 0;
			#print "fail allele : ";
		}
		@line = split(/\s+/,$_);
		if(length($line[3]) > 1 or length($line[4]) > 1){
			$flag = 0;
			#print "fail INDEL : ";
		}
		m/DP=(\d+)/ or die "Syntax error, DP not found\n";
		if ($1 < $minCoverage){
			$flag = 0;
			#print "fail DP : ";
		}
## bcftools call version
	
		m/DP4=\d+,\d+,(\d+),(\d+)/ or die "Syntax error DP4 not found\n";
		if(($1 + $2) < $minAltRds){
			$flag = 0;
		}
		m/AF1*=([0-9\.e\-]+)/ or die "Syntax error, AF not found\n";
		if ($1 == $notFixed){
			$flag = 0;
		#	print "fail AF : ";
		}

## bcftools call verions, these are p-values, use 0.01
		if(m/BQBZ=([0-9e\-\.]*)/){
			if (abs($1) > $bqrs){
				$flag = 0;
#				print "fail BQRS : ";
			}
		}
		if(m/MQBZ=([0-9e\-\.]*)/){
			if (abs($1) > $mqrs){
				$flag = 0;
#				print "fail MQRS : ";
			}
		}
		if(m/RPBZ=([0-9e\-\.]*)/){
			if (abs($1) > $rprs){
				$flag = 0;
#				print "fail RPRS : ";
			}
		}
		if(m/MQ=([0-9\.]+)/){
			if ($1 < $mq){
				$flag = 0;
#				print "fail MQ : ";
			}
		}
		else{
			$flag = 0;
			print "faile no MQ : ";
		}
		if ($flag == 1){
			$cnt++; ## this is a good SNV
		}
	}
	else{
		print "Warning, failed to match the chromosome or scaffold name regular expression for this line\n$_\n";
		$flag = 0;
	}
	if ($flag == 1){
		print OUT "$_\n";
	}
}
close (IN);
close (OUT);

print "Finished filtering $in\nRetained $cnt variable loci\n";
```
Next I extracted the read depth per SNP and individual from the filtered vcf files. I used this information for additional filtering. Specifically, I wanted to remove SNPs with excessively high coverage (possible paralogs) and with a high coefficient of variation for coverage (CV, i.e., high standard deviation in coverage across individuals relative to the mean). I also identified individuals with very high or very low (there were some rather clear high end outliers that seemed worth getting rid of) coverage. This computations were done in `R`, see [CovFilt.R](CovFilt.R), and resulted in vectors of 0s and 1s for SNPs and individuals to drop. For SNPs, I used flagged SNPs to remove with a CV > 1.5 (around the 99.9th percentile) and mean coverage > 3 SDs above the mean (21X per individual). I flagged individuals with mean coverage < 4 (2.5th percentile) or > 40 (a bit above the 90th percentile). I dropped the SNPs first; I will deal with the individuals after the conversion to gl format. I used the following to drop the SNPs, resulting in the morefilter* vcf files.

```bash
#!/bin/bash
#SBATCH --time=240:00:00
#SBATCH --nodes=1
#SBATCH --ntasks=20
#SBATCH --account=gompert-np
#SBATCH --partition=gompert-np
#SBATCH --job-name=filter
#SBATCH --mail-type=FAIL
#SBATCH --mail-user=zach.gompert@usu.edu

cd /scratch/general/nfs1/u6000989/t_chumash_wgs

perl filterSomeMoreL.pl fil*vcf 
```
```perl
#!/usr/bin/perl


# filter vcf files based on coverage


open(IN,"KeepSNPs.txt") or die "failed initial read\n";
while(<IN>){
	chomp;
	push(@keep,$_);
}
close(IN);


foreach $in (@ARGV){
	open (IN, $in) or die "Could not read the infile = $in\n";
	$in =~ m/^([a-zA-Z0-9_]+\.vcf)$/ or die "Failed to match the variant file\n";
	open (OUT, "> morefilter_$1") or die "Could not write the outfile\n";


	while (<IN>){
		chomp;
		if (m/^\#/){ ## header row, always write
			$flag = 1;
		}
		elsif (m/^Sc/){ ## this is a sequence line, you migh need to edit this reg. expr.
			$flag = shift(@keep);
			if ($flag == 1){
				$cnt++; ## this is a good SNV
			}
		}
		else{
			print "Warning, failed to match the chromosome or scaffold name regular expression for this line\n$_\n";
			$flag = 0;
		}
		if ($flag == 1){
			print OUT "$_\n";
		}
	}
	close (IN);
	close (OUT);

	print "Finished filtering $in\nRetained $cnt variable loci\n";
}
```
I then extracted the genotype likelihoods from the filtered vcf files and merged these into a single gl file, called tchum.gl. See [vcf2gl.pl](vcf2gl.pl). This file contained 35,061,459 SNPs.

```bash
## not MAF filter applied, hence 0.0
perl vcf2gl.pl 0.0 morefilter_filtered2x_o_tchum_chrom*vcf
```

# TO DO:

- Split by population or experiment
- EM estimation of allele frequencies
- Empirical Bayes estimate of genotype (mode vs mean?) for experimental population only (for now at least)
- Initial GWA
