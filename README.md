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

I already have the psl file from the alignment of *T. chumash* to the *T. cristinae* melanic, see out_synteny_chumash_melanic.psl. This is what we used for the NEE paper. The genomes I am using beyond this are the *T. chumash*, *T. cristinae* melanic, *T. cristinae* GSH1 and *T. cristinae* GUSH2 (both featured in our recent Science paper). Here is what I have going:

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
I did an initial analysis/investigation of the results (the *T. chumash* ones were actually still running but I the relevant chromosomes were done). The psl files are now in `/uufs/chpc.utah.edu/common/home/gompert-group4/data/timema/hic_genomes/comp_aligns` (but I need to recopy the *T. chumash* ones from scratch once they fully finish). I focused only on chromosome 8 for alignments, see [SynPlotsColorChum.R](SynPlotsColorChum.R). My preliminary plots show that basically have a 1:1 relationship between GSH1 and melanic *T. cristinae* with no compelling evidence of SV. There is an apparent gap (missing scaffold)  between sc1845 and sc128 but this could be simply missing from the melanic assembly not the melanic genome. In contrast, GUSH2 (green haplotype) vs melanic recovers the inverted translocation (see scaffolds 702.1 and 128) (The plot is sort of backwards and upside relative to the plots in the Science paper because the melanic coordinate system is flipped, but the SV looks basically the same, which makes sense given the 1:1 between melanic and striped weware seeing here). Mel-Stripe mostly corresponds with the inverted translocation, though one of the inverted bits isn't fully included. There is no evidence of the deletion from past comparisons between green and melanic, which either means that the coverage and alignment signal resulted from some other form of SV or that the deletion is segregating. For *T. chumash*: (i) the genome is colinear with the striped haplotype rather than the green for the inverted translocation (remeber that stripe matches melanic), (ii) the GWA signal for *T. chumash* overlaps with the translocated inversion, and (iii) it  pecifically overlaps with the colinear bit of the translocation rather than the inverted bit, and it was this colinear bit that was associated with pattern on both refugio and Hwy154.

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

The annotation above resulted in 6151 gene models. I suspect the low number is due to the limited amount of external evidence, especially in terms of RNA sequence data. 

I tried again with additional data from other *Timema* species, specifically 24 RNAseq samples from [Djordjevic et al. 2022](https://www.nature.com/articles/s41437-022-00536-y) (NCBI [PRJNA678950](https://www.ncbi.nlm.nih.gov/sra/?term=PRJNA678950); male and female *T. californicum* hatchling, juveniles and adults) and 232 RNAseq samples from [Djordjevic et al. 2025](https://journals.plos.org/plosgenetics/article?id=10.1371/journal.pgen.1011615) (NCBI [PRJNA1128554](https://www.ncbi.nlm.nih.gov/sra/?term=PRJNA1128554); leg, brain, antennae, gut and gonads from male and female *T. poppense* and *T. douglasi* from multiple developmental stages). I downloaded these `sra-toolkit` to /scratch/general/nfs1/u6000989/rna_expand`, see [grapSRA.pl](grapSRA.pl). 

```perl
#!/usr/bin/perl
#

# ml sra-toolkit

open(IN,"SRAids2.txt") or die;
#open(IN,"SRAids.txt") or die;
while(<IN>){
	chomp;
	system "prefetch $_\n";
}
```

And then extracted the fastq data from the SRA files, see [forkFastqDump.pl](forkFastqDump.pl):

```perl
#!/usr/bin/perl
#
# sra to fastq 
#

# ml sra-toolkit

use Parallel::ForkManager;
my $max = 24;
my $pm = Parallel::ForkManager->new($max);

DIR:
foreach $sra (@ARGV){
	$pm->start and next DIR; ## fork
        system "fasterq-dump ./$sra\n";
	$pm->finish;
}
```
Finally, I then trimmed all of the RNA seqeunces in a consistent manner with `trimmomatic` (version 0.39), see [TrimReadsSimple.sh](TrimReadsSimple.sh).

```bash
#!/bin/sh 
#SBATCH --time=96:00:00
#SBATCH --nodes=1
#SBATCH --ntasks=24
#SBATCH --account=gompert-kp
#SBATCH --partition=gompert-kp
#SBATCH --job-name=trimm
#SBATCH --mail-type=FAIL
#SBATCH --mail-user=zach.gompert@usu.edu


ml trimmomatic/0.39

cd /scratch/general/nfs1/u6000989/rna_expand


# Loop over the input files passed as arguments
for fq1 in "$@"; do
  # Determine the second file by replacing _1.fastq with _2.fastq
  fq2="${fq1/_1.fastq/_2.fastq}"

  # Check if the second file exists
  if [ ! -f "$fq2" ]; then
    echo "Error: $fq2 not found for $fq1. Skipping..."
    continue
  fi

    trimmomatic PE -threads 64 -phred33 "$fq1" "$fq2" \
      "clean_$fq1" "unpaired_$fq1" \
      "clean_$fq2" "unpaired_$fq2" \
      ILLUMINACLIP:Illumina_TruSeq.fa:2:30:8 \
      LEADING:20 TRAILING:20 SLIDINGWINDOW:4:20 MINLEN:35
done
```
I then re-ran the annotation with the added RNA sequence data as shown below.

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
#mkdir ~/augustus
#containerShell
#cp /usr/share/augustus/config ~/augustus
#exit

#cd
#mv ~/augustus/config ~/augustus/config-old
#ml braker/3.0.8
#containerShell
#cp -r /opt/Augustus/config ~/augustus/
#exit

cd /uufs/chpc.utah.edu/common/home/gompert-group4/data/timema/hic_genomes/Annotation/t_chumash

## run braker
braker.pl --genome=/uufs/chpc.utah.edu/common/home/gompert-group4/data/timema/hic_genomes/t_chumash/hic/timema_chumash_29Feb2020_N4ago.Chrom.fasta.masked --prot_seq=../proteins.fasta --rnaseq_sets_ids=clean_SRR13084978,clean_SRR13084979,clean_SRR13084980,clean_SRR13084981,clean_SRR13084982,clean_SRR13084983,clean_SRR13084984,clean_SRR13084985,clean_SRR13084986,clean_SRR13084987,clean_SRR13084988,clean_SRR13084989,clean_SRR13084990,clean_SRR13084991,clean_SRR13084992,clean_SRR13084993,clean_SRR13084994,clean_SRR13084995,clean_SRR13084996,clean_SRR13084997,clean_SRR30782345,clean_SRR30782346,clean_SRR30782347,clean_SRR30782348,clean_SRR30782349,clean_SRR30782350,clean_SRR30782351,clean_SRR30782352,clean_SRR30782353,clean_SRR30782354,clean_SRR30782355,clean_SRR30782356,clean_SRR30782357,clean_SRR30782358,clean_SRR30782359,clean_SRR30782360,clean_SRR30782361,clean_SRR30782362,clean_SRR30782363,clean_SRR30782364,clean_SRR30782365,clean_SRR30782366,clean_SRR30782367,clean_SRR30782368,clean_SRR30782369,clean_SRR30782370,clean_SRR30782371,clean_SRR30782372,clean_SRR30782373,clean_SRR30782374,clean_SRR30782375,clean_SRR30782376,clean_SRR30782377,clean_SRR30782378,clean_SRR30782379,clean_SRR30782380,clean_SRR30782381,clean_SRR30782382,clean_SRR30782383,clean_SRR30782384,clean_SRR30782385,clean_SRR30782386,clean_SRR30782387,clean_SRR30782388,clean_SRR30782389,clean_SRR30782390,clean_SRR30782391,clean_SRR30782392,clean_SRR30782393,clean_SRR30782394,clean_SRR30782395,clean_SRR30782396,clean_SRR30782397,clean_SRR30782398,clean_SRR30782399,clean_SRR30782400,clean_SRR30782401,clean_SRR30782402,clean_SRR30782403,clean_SRR30782404,clean_SRR30782405,clean_SRR30782406,clean_SRR30782407,clean_SRR30782408,clean_SRR30782409,clean_SRR30782410,clean_SRR30782411,clean_SRR30782412,clean_SRR30782413,clean_SRR30782414,clean_SRR30782415,clean_SRR30782416,clean_SRR30782417,clean_SRR30782418,clean_SRR30782419,clean_SRR30782420,clean_SRR30782421,clean_SRR30782422,clean_SRR30782423,clean_SRR30782424,clean_SRR30782425,clean_SRR30782426,clean_SRR30782427,clean_SRR30782428,clean_SRR30782429,clean_SRR30782430,clean_SRR30782431,clean_SRR30782432,clean_SRR30782433,clean_SRR30782434,clean_SRR30782435,clean_SRR30782436,clean_SRR30782437,clean_SRR30782438,clean_SRR30782439,clean_SRR30782440,clean_SRR30782441,clean_SRR30782442,clean_SRR30782443,clean_SRR30782444,clean_SRR30782445,clean_SRR30782446,clean_SRR30782447,clean_SRR30782448,clean_SRR30782449,clean_SRR30782450,clean_SRR30782451,clean_SRR30782452,clean_SRR30782453,clean_SRR30782454,clean_SRR30782455,clean_SRR30782456,clean_SRR30782457,clean_SRR30782458,clean_SRR30782459,clean_SRR30782460,clean_SRR30782461,clean_SRR30782462,clean_SRR30782463,clean_SRR30782464,clean_SRR30782465,clean_SRR30782466,clean_SRR30782467,clean_SRR30782468,clean_SRR30782469,clean_SRR30782470,clean_SRR30782471,clean_SRR30782472,clean_SRR30782473,clean_SRR30782474,clean_SRR30782475,clean_SRR30782476,clean_SRR30782477,clean_SRR30782478,clean_SRR30782479,clean_SRR30782480,clean_SRR30782481,clean_SRR30782482,clean_SRR30782483,clean_SRR30782484,clean_SRR30782485,clean_SRR30782486,clean_SRR30782487,clean_SRR30782488,clean_SRR30782489,clean_SRR30782490,clean_SRR30782491,clean_SRR30782492,clean_SRR30782493,clean_SRR30782494,clean_SRR30782495,clean_SRR30782496,clean_SRR30782497,clean_SRR30782498,clean_SRR30782499,clean_SRR30782500,clean_SRR30782501,clean_SRR30782502,clean_SRR30782503,clean_SRR30782504,clean_SRR30782505,clean_SRR30782506,clean_SRR30782507,clean_SRR30782508,clean_SRR30782509,clean_SRR30782510,clean_SRR30782511,clean_SRR30782512,clean_SRR30782513,clean_SRR30782514,clean_SRR30782515,clean_SRR30782516,clean_SRR30782517,clean_SRR30782518,clean_SRR30782519,clean_SRR30782520,clean_SRR30782521,clean_SRR30782522,clean_SRR30782523,clean_SRR30782524,clean_SRR30782525,clean_SRR30782526,clean_SRR30782527,clean_SRR30782528,clean_SRR30782529,clean_SRR30782530,clean_SRR30782531,clean_SRR30782532,clean_SRR30782533,clean_SRR30782534,clean_SRR30782535,clean_SRR30782536,clean_SRR30782537,clean_SRR30782538,clean_SRR30782539,clean_SRR30782540,clean_SRR30782541,clean_SRR30782542,clean_SRR30782543,clean_SRR30782544,clean_SRR30782545,clean_SRR30782546,clean_SRR30782547,clean_SRR30782548,clean_SRR30782549,clean_SRR30782550,clean_SRR30782551,clean_SRR30782552,clean_SRR30782553,clean_SRR30782554,clean_SRR30782555,clean_SRR30782556,clean_SRR30782557,clean_SRR30782558,clean_SRR30782559,clean_SRR30782560,clean_SRR30782561,clean_SRR30782562,clean_SRR30782563,clean_SRR30782564,clean_SRR30782565,clean_SRR30782566,clean_SRR30782567,clean_SRR30782568,clean_SRR30782569,clean_SRR30782570,clean_SRR30782571,clean_SRR30782572,clean_SRR30782573,clean_SRR30782574,clean_SRR30782575,clean_SRR30782576,clean_tchum_19_0340,clean_tchum_19_0343,clean_tchum_19_0344,clean_tchum_19_0345,clean_tchum_19_0346,clean_tchum_19_0347,clean_tchum_19_0348,clean_tchum_19_0350,clean_tchum_19_0351,clean_tchum_19_0352,clean_tchum_19_0353,clean_tchum_19_0354,clean_tchum_19_0355,clean_tchum_19_0356,clean_tchum_19_0357,clean_tchum_19_0358,clean_tchum_19_0359,clean_tchum_19_0360,clean_tchum_19_0361,clean_tchum_19_0362,clean_tchum_19_0363,clean_tchum_19_0364,clean_tchum_19_0365,clean_tchum_19_0366,clean_tchum_19_0367,clean_tchum_19_0368,clean_tchum_19_0369,clean_tchum_19_0370,clean_tcr135.17_0003_R,clean_tcr137.17_0006_R,clean_tcr139.17_0012_R,clean_tcr140.17_0015_R,clean_tcr141.17_0019_R,clean_tcr142.17_0043_R,clean_tcr143.17_0045_R,clean_tcr144.17_0049_R,clean_tcr145.17_0051_R,clean_tcr146.17_0057_R,clean_tcr148.17_0062_R,clean_tcr149.17_0065_R,clean_tcr150.17_0067_R,clean_tcr151.17_0070_R,clean_tcr152.17_0074_R,clean_tcr173.17_0075_R,clean_tcr174.17_0081_R,clean_tcr175.17_0082_R --rnaseq_sets_dirs=/scratch/general/nfs1/u6000989/rna_expand/ --AUGUSTUS_SCRIPTS_PATH=/usr/share/augustus/scripts --AUGUSTUS_CONFIG_PATH=/uufs/chpc.utah.edu/common/home/u6000989/augustus/config --threads=48 --gff3

## run busco, genome and aa
cd braker

## genome
busco -i /uufs/chpc.utah.edu/common/home/gompert-group4/data/timema/hic_genomes/t_chumash/hic/timema_chumash_29Feb2020_N4ago.Chrom.fasta.masked -m geno -o busco_genome_out -l insecta_odb10

## amino acids
busco -i braker.aa -m prot -o busco_aa_out -l insecta_odb10
```

This gave 20,535 gene models, which is much better and what I will use. The annotation is in: /uufs/chpc.utah.edu/common/home/gompert-group4/data/timema/hic_genomes/Annotation/t_chumash/braker/ (the initial attemp with less RNA data is in brakerv1).

# Looking at candidate genes

In the single SNP GWA, we have 3 peaks of association between ~55 and ~63 Mbp on chromosome 8 (thus spanning a ~8 Mbp region). For the polygenic mapping, the RG signal is at ~59 Mbp and the GB signal is at ~56 and ~59 Mbp. 

I checked out Romain's candidate genes ([Villoutreix et al. 2022](https://royalsocietypublishing.org/doi/full/10.1098/rstb.2020.0508)). We have what looks like two tandem copies of Punch, one shorter than the other. These are at ~34 Mbp on chromosome 8 (about 20 Mbp from the main signal, though there are likely some SNPs with weaker p-value signals closer to these genes). We laso have Chitinase5, which is at ~67.5 Mbp on chromosome 8. None of our annotated genes are good matches for Scarlet (st), though some have a modest similarity. This means either we failed to annotate it (which probably means it was not being expressed in any of the >100 transcriptomes) or it is not present in the T. chumash genome. So, none of his candidates is especially compelling in the current data set.

I also checked our pattern candidates ([Gompert et al. 2025](https://www.science.org/doi/full/10.1126/science.adp3745)). We have a single copy of corazonin, but it is ast ~77 Mbp, so again, somewhat off from the main signal. The only previous candidate withiin the overall range (though not right at 56 or 59 Mbp) is the ecdysteriod kinase.

The main region with the 3 peaks on chromosome 8 (~55 to ~63 Mbp) also contains a cytochrome P450, and those do all kinds of stuff (and have been linked to color, e.g., pigment metaboloism and cuticle tanning; they can even regulate edysteriod metabolism). And right at the far left edge (right where it starts) there is an endonuclease-reverse transcriptase suggestive of a TE insertion (this genes presumably cuts DNA and copies RNA back into the DNA, also interesting int he context of SV). Still, it just seems like an ever growing list of genes and one that isn't really replicating (other than the region, and many genes but not the best candidates, generally coinciding). I was also a bit struck by how big the single SNP GWA signal is: ~8 Mbp.

This brings me to two next steps (with gene expression stuff maybe step 3). First, I want to look at LD across the region and outside of it. I am now wondering about some more minor SV and what to see if there are LD peaks/blocks coinciding with the three peaks. Second, I want to try the polygenic model with only chromosome 8 and see if this alters the SNPs we pick up and the PVE. I want to get a feel for how much to trust them specifically vs the more general three peaks of association. I am also starting to wonder if the causally locus is not a gene at all but rather perhaps some sort of non-coding RNA that regulates something else or something more complicated like that.

# PCA cluster, putative SV, LD and fitness

I finally decided to look for PCA clusters on chromosome 8. We have something that looks like probable SV and that is somewhat associated with color. See (LD.R)[LD.R] and See (LD_GR8.R)[LD_GR8.R]. This led to lots of playing round (need to fill in some details here).

Most important, I looked at the effect of LD on fitness on MM and A/C for the 2019 study population (HFS). This is documented in [Predict.R](Predict.R). I first predicted color (RG and GB) from the gemma results (observed). I then altered LD among top PIP (15 SNPs with PIP > .1 for either color trait) in one of two ways. First, to get minimum LD among these SNPs, I sampled genotypes for individuals independnetly at the 15 loci assuming HWE and LE (that is, binomial sampling based only on allele frequencies). Second, to maximize LD, I used k-means clustering to divide the individuals based on color phenotype (two groups), and then, for each locus, I maximized the allele frequency difference between clusters under the constraint of not altering the allele frequency for both together (basically maximizing genetic differentiation without changing allele frequency for the 15 SNPs). I then sampled genotypes for each locus and individual using the allele frequencies from the relevant cluster/group. Lastly, I caclulated the expected fitness of each individual under the selection gradient models (linear, quadratic and linear correlational selection) from [Nosil et al. 2020](https://www.proquest.com/docview/2473271380?pq-origsite=gscholar&fromopenview=true&sourcetype=Scholarly%20Journals). To do this, I tooks samples from the Bayesian posterior for the higher-level (alpha in the model) selection gradient parameters (as oppossed to block-level parameters). I based on my final inferences on 1000  total samples from the posterior spread over 50 different genotype sampling iterations (20 each). I then computed mean fitness for: (i) observed vs minimum LD = *ran*, (ii) maximum vs observed LD = *max*, and (iii) maximim vs minimum LD = *ext*. The key summaries from [Predict.R](Predict.R) are below. The take home is the LD does more in A/C than MM and we have more credible evidence for an actual beneficial effect of LD in A/C than MM.

```r
relf_ran_MM<-as.vector(obs_mfMM/ran_mfMM)
relf_ran_AC<-as.vector(obs_mfAC/ran_mfAC)
relf_max_MM<-as.vector(max_mfMM/obs_mfMM)
relf_max_AC<-as.vector(max_mfAC/obs_mfAC)
relf_ext_MM<-as.vector(max_mfMM/ran_mfMM)
relf_ext_AC<-as.vector(max_mfAC/ran_mfAC)

## prob > 1
mean(relf_ran_MM > 1)
#[1] 0.812
mean(relf_ran_AC > 1)
#[1] 0.942
mean(relf_max_MM > 1)
#[1] 0.815
mean(relf_max_AC > 1)
#[1] 0.945
mean(relf_ext_MM > 1)
#[1] 0.823
mean(relf_ext_AC > 1)
#[1] 0.984

## median, 95ETPI and 90ETIP
quantile(relf_ran_MM,probs=c(.5,.025,.975,.05,.95))
#      50%      2.5%     97.5%        5%       95% 
#1.0749412 0.9136294 1.2994537 0.9393364 1.2558112 
quantile(relf_ran_AC,probs=c(.5,.025,.975,.05,.95))
#      50%      2.5%     97.5%        5%       95% 
#1.2728799 0.9445857 1.6684976 0.9908804 1.5890200 
quantile(relf_max_MM,probs=c(.5,.025,.975,.05,.95))
#      50%      2.5%     97.5%        5%       95% 
#1.0970248 0.8779211 1.3301346 0.9123359 1.2935863 
quantile(relf_max_AC,probs=c(.5,.025,.975,.05,.95))
#      50%      2.5%     97.5%        5%       95% 
#1.2688999 0.9431728 1.6698536 0.9969875 1.6018711 
quantile(relf_ext_MM,probs=c(.5,.025,.975,.05,.95))
#      50%      2.5%     97.5%        5%       95% 
#1.1897078 0.8087499 1.6783626 0.8716734 1.5878912 
quantile(relf_ext_AC,probs=c(.5,.025,.975,.05,.95))
#     50%     2.5%    97.5%       5%      95% 
#1.595446 1.043744 2.416977 1.148429 2.274322 
```
# Timema cristinae comparative alignments

We are considering adding comparative alignments of all of the phased *T. cristinae* genomes with a focus on chromosomal rearrangements on chromosome 8 associated with color and pattern as part of this project (to think about the evolution of supergenes).

I need to fill in additional details, but, a few quick summaries (before I forget things):

- The new genomes are in `/uufs/chpc.utah.edu/common/home/gompert-group4/data/timema/hic_genomes/edinburgh`

24_0016 VP Green
24_0038 VP Melanic
24_0039 VP Melanic
24_0072 R12 Stripe
24_0073 R12 Green
24_0087 VP Stripe
24_0175 FH Stripe
24_0176 FH Stripe
  
- I ran repeat masking on them and did a big series of preliminary alignments, see [SynPlotTcrEd.R](SynPlotTcrEd.R)

- I then extracted chromosome 8 from all of the phased, *T. cristinae* genomes, see [extractCh8.pl](extractCh8.pl) and `/uufs/chpc.utah.edu/common/home/gompert-group4/data/timema/hic_genomes/comp_aligns/chr8haplotypes`

- I used SibeliaZ to do an initial alignment of this (it is fast but not as good at deeper divergence)

```bash
#!/bin/bash 
#SBATCH --time=96:00:00
#SBATCH --nodes=1
#SBATCH --ntasks=96
#SBATCH --mem=724000
#SBATCH --account=gompert
#SBATCH --qos=gompert-grn
#SBATCH --partition=gompert-grn
#SBATCH --job-name=sibelia
#SBATCH --mail-type=FAIL
#SBATCH --mail-user=zach.gompert@usu.edu

## sibelia comparative alignment for ch8

cd /uufs/chpc.utah.edu/common/home/gompert-group4/data/timema/hic_genomes/comp_aligns/chr8haplotypes

## see https://github.com/medvedevgroup/SibeliaZ
## for options
## -a = 24 genomes * 4 * 2
~/source/SibeliaZ/build/bin/bin/sibeliaz -k 25 -a 192 -b 200 -m 50 -t 180 ch8_t_cris_e_240016h1.fasta ch8_t_cris_e_240072h1.fasta ch8_t_cris_e_240175h1.fasta ch8_t_cris_h_gus1.fasta ch8_t_cris_e_240016h2.fasta ch8_t_cris_e_240072h2.fasta ch8_t_cris_e_240175h2.fasta ch8_t_cris_h_gus2.fasta ch8_t_cris_e_240038h1.fasta ch8_t_cris_e_240073h1.fasta ch8_t_cris_e_240176h1.fasta ch8_t_cris_r_gs1.fasta ch8_t_cris_e_240038h2.fasta ch8_t_cris_e_240073h2.fasta ch8_t_cris_e_240176h2.fasta ch8_t_cris_r_gs2.fasta ch8_t_cris_e_240039h1.fasta ch8_t_cris_e_240087h1.fasta ch8_t_cris_h_gs1.fasta ch8_t_cris_r_gus1.fasta ch8_t_cris_e_240039h2.fasta ch8_t_cris_e_240087h2.fasta ch8_t_cris_h_gs2.fasta ch8_t_cris_r_gus2.fasta
```
- The results are in `/uufs/chpc.utah.edu/common/home/gompert-group4/data/timema/hic_genomes/comp_aligns/chr8haplotypes/sibeliaz_out`

- I converted the aligment file to mfa format using [maftools](https://github.com/dentearl/mafTools/tree/master?tab=readme-ov-file)

```bash
 ~/source/mafTools/bin/mafToFastaStitcher --maf sibeliaz_out/alignment.maf --seqs ch8_combined.fasta --breakpointPenalty 5 --interstitialSequence 20  --outMfa output.mfa
```

- Next I used [trimAl](https://trimal.readthedocs.io/en/latest/usage.html) to remove any positions from the alignment with gaps.

```bash
## auto, not used
## trimal -in output.mfa -out clean.mfa -gappyout
## remove positions with any gaps
trimal -in output.mfa -out clean.mfa -gt 1.0
```

- I used castertree topology with caster

```bash
~/../gompert-group5/projects/LycAdmix/Caster/ASTER-Linux/bin/caster-site -i clean.mfa -o cout_tree --thread 24
```

- Finally, I added ML branch lengths to the tree in R, see [AddBlen.R](AddBlen.R).

```R
## adds ml branch lengths to the tree
library(phangorn)
library(ape)

packageVersion("phangorn")
#[1] ‘2.12.1’
packageVersion("ape")
#[1] ‘5.8’

## read tree, just topology, from caster
tre<-read.tree("topo_cout_tree")

## read the alignment and convert to phyDat object
aln<-read.dna("clean.mfa",format="fasta")
alnPh<-phyDat(aln,type="DNA")

## add arbitrary branch lengths as a starting point
treb<-compute.brlen(tre,1)

## compute the initial likelihood and then optimize branch lenghts under JC
## given the topoloty
fit<-pml(tree=treb,data=alnPh)
fit_optimized <- optim.pml(fit, model = "JC", optEdge = TRUE, optGamma = FALSE, optInv = FALSE)

## I might want to increase these slightly so they are overestimates
fit_optimized$tree$edge.length<-fit_optimized$tree$edge.length * 1.5

## write the tree
write.tree(phy=fit_optimized$tree,file="topoB_cout_tree")

## root on refugio and write again
refugio<-fit_optimized$tree$tip.label[9:16]
trr<-root(fit_optimized$tree,outgroup=refugio)
write.tree(phy=trr,file="topoRoot_cout_tree")
```

- Lastly, I used this to make the Cactus infile, I had to manually add the ancestors and fix the Refugio root

```
(((((((t_cris_e_240087h2:0.4426538755,t_cris_h_gs2:0.4582198529):0.07364615261,t_cris_e_240087h1:0.4635711029)AVPS:0.1559547772,(((t_cris_e_240016h2:0.3290234023,t_cris_h_gs1:0.2991242597):0.04056500028,t_cris_e_240039h1:0.4349549599):0.03885041141,(t_cris_e_240038h1:0.292011903,t_cris_e_240038h2:0.3390838267):0.06219885817)AVPM:0.03646449823):0.06323898609,t_cris_e_240039h2:0.4271548927):0.02428566702,(((t_cris_e_240176h1:0.304560471,t_cris_e_240175h1:0.3542393032):0.06988151601,t_cris_e_240175h2:0.2887467482):0.09896247756,t_cris_e_240176h2:0.3692333657)AFHS:0.05223943831):0.032424558,((t_cris_h_gus1:0.2903114095,t_cris_h_gus2:0.4758871493):0.1156325651,t_cris_e_240016h1:0.4679979238)AVPG:0.05056619986)AH154:0.05505660008,(((((t_cris_r_gus2:0.3435745311,t_cris_r_gus1:0.3564212078):0.1041618822,t_cris_e_240073h1:0.3693590572):0.04537837434,t_cris_e_240073h2:0.415594649)ARG:0.07449178421,((t_cris_r_gs2:0.346746946,t_cris_e_240072h1:0.3992755297):0.06593359072,t_cris_r_gs1:0.2781758046)ARS:0.1200039446):0.1166655065,t_cris_e_240072h2)AREF:0.4529215631);

t_cris_e_240016h1 /uufs/chpc.utah.edu/common/home/gompert-group4/data/timema/hic_genomes/edingburgh/24_0016/Hap1Chr.fasta.masked
t_cris_e_240016h2 /uufs/chpc.utah.edu/common/home/gompert-group4/data/timema/hic_genomes/edingburgh/24_0016/Hap2Chr.fasta.masked
t_cris_e_240038h1 /uufs/chpc.utah.edu/common/home/gompert-group4/data/timema/hic_genomes/edingburgh/24_0038/Hap1Chr.fasta.masked
t_cris_e_240038h2 /uufs/chpc.utah.edu/common/home/gompert-group4/data/timema/hic_genomes/edingburgh/24_0038/Hap2Chr.fasta.masked
t_cris_e_240039h1 /uufs/chpc.utah.edu/common/home/gompert-group4/data/timema/hic_genomes/edingburgh/24_0039/Hap1Chr.fasta.masked
t_cris_e_240039h2 /uufs/chpc.utah.edu/common/home/gompert-group4/data/timema/hic_genomes/edingburgh/24_0039/Hap2Chr.fasta.masked
t_cris_e_240072h1 /uufs/chpc.utah.edu/common/home/gompert-group4/data/timema/hic_genomes/edingburgh/24_0072/Hap1Chr.fasta.masked
t_cris_e_240072h2 /uufs/chpc.utah.edu/common/home/gompert-group4/data/timema/hic_genomes/edingburgh/24_0072/Hap2Chr.fasta.masked
t_cris_e_240073h1 /uufs/chpc.utah.edu/common/home/gompert-group4/data/timema/hic_genomes/edingburgh/24_0073/Hap1Chr.fasta.masked
t_cris_e_240073h2 /uufs/chpc.utah.edu/common/home/gompert-group4/data/timema/hic_genomes/edingburgh/24_0073/Hap2Chr.fasta.masked
t_cris_e_240087h1 /uufs/chpc.utah.edu/common/home/gompert-group4/data/timema/hic_genomes/edingburgh/24_0087/Hap1Chr.fasta.masked
t_cris_e_240087h2 /uufs/chpc.utah.edu/common/home/gompert-group4/data/timema/hic_genomes/edingburgh/24_0087/Hap2Chr.fasta.masked
t_cris_e_240175h1 /uufs/chpc.utah.edu/common/home/gompert-group4/data/timema/hic_genomes/edingburgh/24_0175/Hap1Chr.fasta.masked
t_cris_e_240175h2 /uufs/chpc.utah.edu/common/home/gompert-group4/data/timema/hic_genomes/edingburgh/24_0175/Hap2Chr.fasta.masked
t_cris_e_240176h1 /uufs/chpc.utah.edu/common/home/gompert-group4/data/timema/hic_genomes/edingburgh/24_0176/Hap1Chr.fasta.masked
t_cris_e_240176h2 /uufs/chpc.utah.edu/common/home/gompert-group4/data/timema/hic_genomes/edingburgh/24_0176/Hap2Chr.fasta.masked
t_cris_h_gs1 /uufs/chpc.utah.edu/common/home/gompert-group4/data/timema/hic_genomes/t_crist_gs_hap_cen4119/HiRise/Hap1/chroms_final_assembly.fasta.masked
t_cris_h_gs2 /uufs/chpc.utah.edu/common/home/gompert-group4/data/timema/hic_genomes/t_crist_gs_hap_cen4119/HiRise/Hap2/chroms_final_assembly.fasta.masked
t_cris_h_gus1 /uufs/chpc.utah.edu/common/home/gompert-group4/data/timema/hic_genomes/t_crist_gus_hap_cen4280/HiRise/Hap1/chroms_final_assembly.fasta.masked
t_cris_h_gus2 /uufs/chpc.utah.edu/common/home/gompert-group4/data/timema/hic_genomes/t_crist_gus_hap_cen4280/HiRise/Hap2/chroms_final_assembly.fasta.masked
t_cris_r_gs1 /uufs/chpc.utah.edu/common/home/gompert-group4/data/timema/hic_genomes/t_crist_refug_stripe/HiRise/hap1/chroms_final_assembly.fasta.masked
t_cris_r_gs2 /uufs/chpc.utah.edu/common/home/gompert-group4/data/timema/hic_genomes/t_crist_refug_stripe/HiRise/hap2/chroms_final_assembly.fasta.masked
t_cris_r_gus1 /uufs/chpc.utah.edu/common/home/gompert-group4/data/timema/hic_genomes/t_crist_refug_green/HiRise/hap1/chroms_final_assembly.fasta.masked
t_cris_r_gus2 /uufs/chpc.utah.edu/common/home/gompert-group4/data/timema/hic_genomes/t_crist_refug_green/HiRise/hap2/chroms_final_assembly.fasta.masked
```
```bash
#!/bin/bash 
#SBATCH --time=168:00:00
#SBATCH --nodes=1
#SBATCH --ntasks=96
#SBATCH --mem=524000
#SBATCH --account=gompert
#SBATCH --qos=gompert-grn
#SBATCH --partition=gompert-grn
#SBATCH --job-name=cactus
#SBATCH --mail-type=FAIL
#SBATCH --mail-user=zach.gompert@usu.edu

module load cactus

cd /scratch/general/nfs1/u6000989/cactus

cactus jobStore_Prog /uufs/chpc.utah.edu/common/home/gompert-group4/data/timema/hic_genomes/comp_aligns/chr8haplotypes/ch8_cactus.txt cactusTcrAll8.hal --maxCores 72
```
