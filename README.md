# ChumashWGSmapping
Genomic analyses of color and survival in a selection experiment for *Timema chumash* based on whole genome sequence data

# Data

The DNA sequence data for this project are in two places. The first two batches are in `/uufs/chpc.utah.edu/common/home/gompert-group4/data/t_chumash_wgs/`: batch 1 = 192 samples, batch 2 = 288 samples (some of the latter are not part of the mapping and selection experiment). See [Nosil et al. 2020](https://www.proquest.com/docview/2473271380?pq-origsite=gscholar&fromopenview=true&sourcetype=Scholarly%20Journals) ([Nosil2020.pdf](https://github.com/user-attachments/files/17368867/Nosil2020.pdf)) for our analyses of GBS data from this same experiment. The data from the 3rd batch are in `/uufs/chpc.utah.edu/common/home/gompert-group5/data/TchumJan2025/`. Batch 3 includes 541 additioanl *T. chumash*

**batch 1:** These samples were sequenced by GeT (Genotoul) and delivered July 2024. The IDs are at the start of the file names (Tchum-XXX, where XXX denotes the digits after 19_ in the actual file names). There are two files each of forward and reverse reads for each stick insect. 

**batch 2:** These samples were sequenced by GeT (Genotoul) and delivered October 2024. The IDs are again at the start of the file name but here are either 19_XXX (the actual ID we use) or GR806-XX for the non-experimental stick insects (the GR806-XX samples were not in a spreadsheet but come from a sample in ethanol from host MM).

**batch 3:** These samples were sequenced by GeT (Genotoul) and delivered January 2025. The IDs are again at the start of the file name but here are either 15_XXX. These are all non-experimental *T. chumash* from natural populations. The population information is in `/uufs/chpc.utah.edu/common/home/gompert-group5/data/TchumJan2025/sample_sp_loc_host.dsv'. The GR806-XX samples from batch 2 go with these individuals (i.e., they add to this data set).

The experimental data (color and survival) are in the 2019_Tchumash* files in `/uufs/chpc.utah.edu/common/home/gompert-group4/data/t_chumash_wgs/`.

I am using our 2020 *T. chumash* genome assembly. This genome is described in [Nosil et al. 2020](https://www.proquest.com/docview/2473271380?pq-origsite=gscholar&fromopenview=true&sourcetype=Scholarly%20Journals). See [timema_chumash_29Feb2020_N4ago.report.pdf](https://github.com/user-attachments/files/17368853/timema_chumash_29Feb2020_N4ago.report.pdf) for a report from Dovetail. Past comparative alignments suggest general chromosome-level synteny with *T. crstinae* except that *T. cristinae* chromsomes 1, 3, 4 and 9 appear to be fused in *T. chumash*. Patrik suggested this is consistent with cytology (add source). Also note (from the report) that the genome size is notably larger than our other *Timema* genomes (2.4 Gbps) and the final assembly has a non-trivial number of duplicated BUSCOs (29). This could reflect some difficulities in the assembly (especially if the two genome copies differed enough) that might be cleaned up with phased genomes later. However, chromosome 8 appears to be largely syntenic (though not colinear) between *T. cristinae* and *T. chumash* and is the main focus of our interest ([see the TimemaFusion repository for past work on this](https://github.com/zgompert/TimemaFusion)).

We have color data for the 2019 and 2015 (GR806) samples; we have survival data for 2019 only (see [2019_Tchumash_transplant_phenotype.csv](2019_Tchumash_transplant_phenotype.csv), [2019_Tchumash_transplant_table.csv](2019_Tchumash_transplant_table.csv) and [2015_Tchumash_color.csv](2015_Tchumash_color.csv)). I had the 2019 trait data on the cluster and obtained the 2015 data from Dryad. Both are now also in `/uufs/chpc.utah.edu/common/home/gompert-group5/projects/t_chum_mapping/gemma`.

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

# Genotype and allele frequency estimation

We had 35,061,459 in the initial genotype likelihood file, tchum.gl. I split this initially into files by population (locality and host). I then went back and did a second split to have a combiend GR8.06 (not by Q and MM) as this will be useful for mapping (I still might want by host for allele frequencies, etc.). As part of this process, I removed individuals with low coverage (based on depth.txt, depth extracted from the vcf files) and flagged SNPs with low coverage to remove later. This was done with [CovFilt.R](CovFilt.R). This wrote the files KeepInds.txt and KeepSNPs.txt. Low coverage individuals were split into a NA population (for the trash). 

```bash
#!/bin/bash
#SBATCH --time=240:00:00
#SBATCH --nodes=1
#SBATCH --ntasks=12
#SBATCH --account=gompert-kp
#SBATCH --partition=gompert-kp
#SBATCH --job-name=glsplit
#SBATCH --mail-type=FAIL
#SBATCH --mail-user=zach.gompert@usu.edu

cd /scratch/general/nfs1/u6000989/t_chumash_wgs

perl splitPops.pl tchum.gl 
```
Which runs [splitPops.pl](splitPops.pl).

The GR8.06 file was made with [splitPops2.pl](splitPops2.pl).

Next, I used my expectation-maximization program, `estpEM` (version 0.1) (from [Soria-Carrasco et al 204](https://www.science.org/doi/full/10.1126/science.1252136)), to estimate population allele frequencies.  

```bash
estpEM -i EXP_tchum.gl -o p_EXP_tchum.txt -e 0.001 -m 50 -h 2 
estpEM -i tchum.gl -o p_tchum.txt -e 0.001 -m 50 -h 2 
estpEM -i GR8.06_MM_tchum.gl -o p_GR8.06_MM_tchum.txt -e 0.001 -m 50 -h 2 
estpEM -i GR8.06_Q_tchum.gl -o p_GR8.06_Q_tchum.txt -e 0.001 -m 50 -h 2
estpEM -i GR8.06_tchum.gl -o p_GR8.06_tchum.txt -e 0.001 -m 50 -h 2 
```
These files include allele frquency estimates for 35,061,459 loci, with 429 individuals from the 2019 experiment (Horse Flats) and 555 from the 2015 GR8.06 population.

I then obtained Bayesian estiamtes of genotypes using the allele frequency priors (under assuming HW genotype frequencies as prior expectations). I did this using both the posterior mode and mean. I am processing the results from the mode first, but want to see how much this matters. I used new C programs, based on my older perl script, for the empirical Bayes genotype esimates, see [gl2genest.c](gl2genest.c) and [gl2genestMax.c](gl2genestMax.c).

```bash
#!/bin/bash
#SBATCH --time=48:00:00
#SBATCH --nodes=1
#SBATCH --ntasks=12
#SBATCH --account=gompert
#SBATCH --partition=kingspeak
#SBATCH --job-name=gl2g
#SBATCH --mail-type=FAIL
#SBATCH --mail-user=zach.gompert@usu.edu

cd /scratch/general/nfs1/u6000989/t_chumash_wgs

## posterior mode
./gl2genestMax pp_EXP_tchum.txt EXP_tchum.gl
./gl2genestMax pp_GR8.06_tchum.txt GR8.06_tchum.gl
## posterior mean
./gl2genest pp_EXP_tchum.txt EXP_tchum.gl
./gl2genest pp_GR8.06_tchum.txt GR8.06_tchum.gl
```
The output files are cpntest_EXP_tchum.txt for the posterior mean and mlpntest_EXP_tchum.txt for the posterior mode. These are for the 2019 experiment and contain 429 individuals (columns) and 35,061,459 SNPs (rows). These files, along with the gl files, are in `/uufs/chpc.utah.edu/common/home/gompert-group5/projects/t_chum_mapping/gendat/`.


# Genome-wide association mapping of color (and survival?)

Next, I prepared the genetic data for GWA mapping with `gemma`. I did this first for the 2019 data set. I first identified common SNPs within that data set, that is SNPs with a minor allele frequency > 0.01. See [getCommon.R](getCommon.R). I also extracted the scaffold, position and allele information from the vcf files, which I wrote to a file = FullSNPTable.txt. This and the set of common SNPs (KeepSnps.txt) are in `/uufs/chpc.utah.edu/common/home/gompert-group5/projects/t_chum_mapping/gendat`. I then ran [FormatGeno.pl](FormatGeno.pl) to create the geno file from the model genotype esimate file.

```bash
#!/bin/bash
#SBATCH --time=48:00:00
#SBATCH --nodes=1
#SBATCH --ntasks=12
#SBATCH --account=gompert
#SBATCH --partition=notchpeak
#SBATCH --job-name=vcf2gl
#SBATCH --mail-type=FAIL
#SBATCH --mail-user=zach.gompert@usu.edu

cd /scratch/general/nfs1/u6000989/t_chumash_wgs

perl FormatGeno.pl KeepSnps.txt FullSNPTable.txt mlpntest_EXP_tchum.txt
```

```perl
#!/usr/bin/perl
#
# takes a file with a vector of SNPs to keep
# a snp info file and a genotype fileand generates gemma input
# example: CommonSnps.txt FullSNPTable.txt mlpntest_EXP_tchum.txt

$kf = shift(@ARGV); ## keep file
$sf = shift(@ARGV); ## SNP file
$gf = shift(@ARGV); ## genotype file

## get set of SNPs to keep
open(IN, $kf) or die "failed to read $kf\n";
while(<IN>){
	chomp;
	push(@keep,$_);
}
close(IN);

## get SNP information
open(IN, $sf) or die;
while(<IN>){
	chomp;
	s/^morefilter_filtered2x_o_tchum_chrom(\d+)\S+\s+/\1:/;
	s/\s+/ /g;
	push(@snps,$_);
}
close(IN);

print "Formatting genotype file now\n";

## subset and format genotype file
open(IN, $gf);
$out = $gf;
$out =~ s/txt/geno/ or die "failed sub: $out\n";
open(OUT, "> $out") or die;
while(<IN>){
	chomp;
	$k = shift(@keep);
	$s = shift(@snps);
	if($k == 1){
		$o = "$s $_";
		$o =~ s/ /, /g;
		print OUT "$o\n";
	}
}
close(IN);
close(OUT);
```
Next, working in `/uufs/chpc.utah.edu/common/home/gompert-group5/projects/t_chum_mapping/gemma`, I preprated the phenotypic data for the 2019 experimental popualtion, this includes RG, GB and (binary) survival. Right now survival is across both host treatments (I will split this by host later as that makes the most sense). See [FormatPheno.R](FormatPheno.R).

I also got data on observed heterzygosity, deviations between observed and expected heterozygosity and allele frequencies to filter out potentially problematic SNPs. See [GetHetPdata.R](GetHetPdata.R). I used this information to create a subset geno file that keeps only SNPs with deviations between observed and expected heterozygosity < .2. See [SubSetHetGeno.R](SubSetHetGeno.R). I ran initial LMM analyses on color and surival on the geno file before subsetting, though I then did some subsetting in plots after the fact. 

For the LMM.

```bash
#!/bin/bash
#SBATCH --time=96:00:00
#SBATCH --nodes=1
#SBATCH --ntasks=12
#SBATCH --account=wolf-kp
#SBATCH --partition=wolf-kp
#SBATCH --job-name=gemma
#SBATCH --mail-type=FAIL
#SBATCH --mail-user=zach.gompert@usu.edu

module load gemma
# 0.95a

cd /uufs/chpc.utah.edu/common/home/gompert-group5/projects/t_chum_mapping/gemma

## RG 
gemma -g  mlpntest_EXP_tchum.geno -p ph_tchum2019.txt -gk 1 -o o_RG -maf 0
gemma -g  mlpntest_EXP_tchum.geno -p ph_tchum2019.txt -k output/o_RG.cXX.txt -lmm 4 -n 1 -o o_RG -maf 0

## GB 
gemma -g  mlpntest_EXP_tchum.geno -p ph_tchum2019.txt -gk 1 -o o_GB -maf 0
gemma -g  mlpntest_EXP_tchum.geno -p ph_tchum2019.txt -k output/o_GB.cXX.txt -lmm 4 -n 2 -o o_GB -maf 0

## Surv 
#gemma -g  mlpntest_EXP_tchum.geno -p ph_tchum2019.txt -gk 1 -o o_SURV -maf 0
gemma -g  mlpntest_EXP_tchum.geno -p ph_tchum2019.txt -k output/o_SURV.cXX.txt -lmm 4 -n 3 -o o_SURV -maf 0
```
This was summarized with [summarizeGwa.R](summarizeGwa.R).

Next, I ran the BSLMM with the subset version. I ran 60 chains total.

```bash
#!/bin/bash
#SBATCH --time=240:00:00
#SBATCH --nodes=1
#SBATCH --ntasks=12
#SBATCH --account=wolf-kp
#SBATCH --partition=wolf-kp
#SBATCH --job-name=gemma
#SBATCH --mail-type=FAIL
#SBATCH --mail-user=zach.gompert@usu.edu

module load gemma
# 0.95a

cd /uufs/chpc.utah.edu/common/home/gompert-group5/projects/t_chum_mapping/gemma

## RG 
gemma -g  hsub_mlpntest_EXP_tchum.geno -p ph_tchum2019.txt -gk 1 -o o_poly_RG -maf 0
gemma -g  hsub_mlpntest_EXP_tchum.geno -p ph_tchum2019.txt -k output/o_poly_RG.cXX.txt -lmm 4 -n 1 -o o_poly_RG -maf 0

## GB 
gemma -g  hsub_mlpntest_EXP_tchum.geno -p ph_tchum2019.txt -gk 1 -o o_poly_GB -maf 0
gemma -g  hsub_mlpntest_EXP_tchum.geno -p ph_tchum2019.txt -k output/o_poly_GB.cXX.txt -lmm 4 -n 2 -o o_poly_GB -maf 0

## Surv 
gemma -g  hsub_mlpntest_EXP_tchum.geno -p ph_tchum2019.txt -gk 1 -o o_poly_SURV -maf 0
gemma -g  hsub_mlpntest_EXP_tchum.geno -p ph_tchum2019.txt -k output/o_poly_SURV.cXX.txt -lmm 4 -n 3 -o o_poly_SURV -maf 0

## now polygenic BSLMM
perl forkRunGemma.pl
```
The perl script is:
```perl
#!/usr/bin/perl
#
# fit gemma BSLMM for T. chumash color 
#

use Parallel::ForkManager;
my $max = 40;
my $pm = Parallel::ForkManager->new($max);

$g = "hsub_mlpntest_EXP_tchum.geno";
$p = "ph_tchum2019.txt";


foreach $ph (1..2){ 
	foreach $ch (0..59){
		sleep 2;
		$pm->start and next;
		$o = "o_poly_bslmm_ph$ph"."_ch$ch";
		if($ph == 1){
	    		system "gemma -g $g -p $p -bslmm 1 -n $ph -o $o -maf 0 -w 200000 -s 1000000 -k output/o_poly_RG.cXX.txt\n";
		} elsif($ph == 2){
	    		system "gemma -g $g -p $p -bslmm 1 -n $ph -o $o -maf 0 -w 200000 -s 1000000 -k output/o_poly_GB.cXX.txt\n";
		}
		$pm->finish;
	}
}
$pm->wait_all_children;
```

I summarized the results from the chains, see [calpost.pl](calpost.pl) and [grabPips.pl](grabPips.pl), and then summarized the mapping analysis in R [summarizeGwaPoly.R](summarizeGwaPoly.R).

I next repeated all of the above (but with subsetting first) for GR8.06. See [FormatPhenoGR8.R](FormatPhenoGR8.R), [GetHetPdataGR8.R](GetHetPdataGR8.R), [SubSetHetGenoGR8.R](SubSetHetGenoGR8.R), [summarizeGwaGR.R](summarizeGwaGR.R), [grabPipsGR8.pl](grabPipsGR8.pl) and [summarizeGwaPolyGR.R](summarizeGwaPolyGR.R). 

I put all of the figures we have so far on [google drive folder](https://drive.google.com/drive/folders/1eIxAUx4Up3Opkr4DiMcKhPcZjayXjwke?usp=sharing). The file names denote the population (GR8 for GR8.06 and 2019 for Horse Flats), the analysis (SSNP = single SNP LMM vs poply = BSLMM), the trait (GB or GR) and whether the plot is for the whole genome (nothing extra in the name) or just chromosome 8 (Ch8). All combinations of these are included for 16 plots (files total).

We see fairly similar results in both populations. PVEs (with CIs) are: Horse Flats, RG = 0.75 (0.57-0.92), GB = 0.63 (0.47-0.82); GR8.06, RG = 0.72 (0.60-0.91), GB = 0.60 (0.48-0.76). In both cases the color traits are highly heritable, with slightly higher heritability for RG than GB and similar values in each population.

In all cases we see a clear signal on chromosome 8. There is generally more noise off of 8 (sometimes quite a bit) with the single SNP analyses than the polygenic analyses. In the polygenic analyses, we get a modest number of SNPs with non-trivial PIPs (more details below). For GB in both populations, the signal in the polygenic analyses is very much restricted to chromosome 8. For RG, there is a clear signal on 8 in the polygenic models, but also something going on with chromosome 1 (both populations) with other more modest signals elsewhere.

The overall (total) number of sparse effects (QTN) from gemma is modest and exhibits uncertainty (of course) but not to a crazy extent: Horse Flats, RG = 18 (8-39), GB = 15 (5-53); GR8.06, RG = 10 (7-18), GB = 12 (3-35). Interestingly, there is some evidence that sparse effects matter more for GR8.06 than Horse Flats, though it isn't whopping, based on PGE (proportion of PVE due to sparse effects): Horse Flats, RG = 0.79 (0.60-0.96), GB = 0.77 (0.54-0.97); GR8.06, RG = 0.88 (0.70-0.99), GB = 0.90 (0.70-0.99). This combined with the slightly lower estimates of number of QTN point towards a marginally more concentrated (less polygenic) genetic basis for color at GR8.06 than Horse Flats. The SNPs with the highest PIPs are not the same for RG and GB, but, at least on chromosome 8, they are in the same general region (don't know about gene level yet), at least for the ones with the strongest signal. It also looks like, especially when thinking about RG and GB, we get a few distinct, but nearby peaks. This is all generally consistent with the NEE paper.

Based on the PIPs, we have 3.2 and 3.6 QTN for RG and GB (respectively) on Ch8 in Horse Flats and 3.3 and 3.4 QTN for RG and GB (respectively) on Ch8 in GR8.06. If we only consider PIPs > 0.01 as contributing to these estimates (which gets rid of some of the noise) this drops to 2.2 and 2.2 QTN for RG and GB (respectively) on Ch8 in Horse Flats and 2.9 and 2.5 QTN for RG and GB (respectively) on Ch8 in GR8.06. This contains almost all such SNPS for GB, but more like half of such SNPs for RG (where there is more signal off 8 in the polygenic models).

# Comparative genome alignments, *T. chumash* with old melanic and phased *T. cristinae*

Here, I am mosly interested in cross-referencing the old color signal (from melanic deletion thinking) with the new *T. cristiane* genomes and the new *T. chumash* mapping results, which are based on the *T. chumash* genome coordinate system.

I already have the psl file from the alignment of *T. chumash* to the *T. cristinae* melanic, see out_synteny_chumash_melanic.psl. This is what we used for the NEE paper. The genomes I am using beyond this are the *T. chumash*, T. cristinae* melanic, *T. cristinae* GSH1 and *T. cristinae* GUSH2 (both featured in our recent Science paper). Here is what I have going:

```bash
#!/bin/bash 
#SBATCH --time=240:00:00
#SBATCH --nodes=1
#SBATCH --ntasks=18
#SBATCH --account=gompert-np
#SBATCH --partition=gompert-np
#SBATCH --job-name=cactus
#SBATCH --mail-type=FAIL
#SBATCH --mail-user=zach.gompert@usu.edu

module load cactus

cd /scratch/general/nfs1/u6000989/cactus

perl CactusForkColor.pl
```
Which runs:

```perl
#!/usr/bin/perl
#
# cactus batch run 
#


use Parallel::ForkManager;
my $max = 4;
my $pm = Parallel::ForkManager->new($max);


@comps = ("cactusColor_TcrGSH1_Tchum.txt","cactusColor_TcrGUSH2_Tchum.txt","cactusColor_TcrGSH1_TcrM.txt","cactusColor_TcrGUSH2_TcrM.tx");

@ids1 = ("t_cris_h_gs1","t_cris_h_gus2","t_cris_h_gs1","t_cris_h_gus2");
@ids2 = ("t_chum","t_chum","t_cris_m","t_cris_m");

FILES:
foreach $file (@comps){
	$id1 = shift(@ids1);
	$id2 = shift(@ids2);
	$pm->start and next FILES; ## fork
	$file =~ m/cactusColor_([a-zA-Z0-9]+)_([a-zA-Z]+)/ or die "failed here $file\n";
	$base = "$1_$2";
	system "cactus jobStore_$base /uufs/chpc.utah.edu/common/home/gompert-group4/data/timema/hic_genomes/comp_aligns/$file cactus$base.hal --maxCores 24\n";
	system "~/source/hal/bin/halSynteny --queryGenome $id2 --targetGenome $id1 cactus$base.hal out_synteny_$base.psl\n";
	
	$pm->finish;
}

$pm->wait_all_children;
```

## *T. chumash* genome annotation with BRAKER3

```bash
#!/bin/bash 
#SBATCH --time=240:00:00
#SBATCH --nodes=1
#SBATCH --ntasks=24
#SBATCH --account=gompert-np
#SBATCH --partition=gompert-np
#SBATCH --job-name=braker
#SBATCH --mail-type=FAIL
#SBATCH --mail-user=zach.gompert@usu.edu

source ~/.bashrc

ml braker/3.0.8
ml busco

cd /uufs/chpc.utah.edu/common/home/gompert-group4/data/timema/hic_genomes/Annotation/t_chumash

## run braker
braker.pl --genome=/uufs/chpc.utah.edu/common/home/gompert-group4/data/timema/hic_genomes/t_chumash/hic/timema_chumash_29Feb2020_N4ago.Chrom.fasta.masked --prot_seq=../proteins.fasta --rnaseq_sets_ids=tcr135.17_0003_R,tcr137.17_0006_R,tcr139.17_0012_R,tcr140.17_0015_R,tcr141.17_0019_R,tcr142.17_0043_R,tcr143.17_0045_R,tcr144.17_0049_R,tcr145.17_0051_R,tcr146.17_0057_R,tcr148.17_0062_R,tcr149.17_0065_R,tcr150.17_0067_R,tcr151.17_0070_R,tcr152.17_0074_R,tcr173.17_0075_R,tcr174.17_0081_R,tcr175.17_0082_R --rnaseq_sets_dirs=/scratch/general/nfs1/u6000989/rna/ --AUGUSTUS_SCRIPTS_PATH=/usr/share/augustus/scripts --AUGUSTUS_CONFIG_PATH=/uufs/chpc.utah.edu/common/home/u6000989/augustus/config --threads=48 --gff3

## run busco, genome and aa
cd braker

## genome
busco -i /uufs/chpc.utah.edu/common/home/gompert-group4/data/timema/hic_genomes/t_chumash/hic/timema_chumash_29Feb2020_N4ago.Chrom.fasta.masked -m geno -o busco_genome_out -l insecta_odb10

## amino acids
busco -i braker.aa -m prot -o busco_aa_out -l insecta_odb10
```
