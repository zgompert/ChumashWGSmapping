# ChumashWGSmapping
Genomic analyses of color and survival in a selection experiment for Timema chumash based on whole genome sequence data

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




