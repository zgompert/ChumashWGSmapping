# ChumashWGSmapping
Genomic analyses of color and survival in a selection experiment for Timema chumash based on whole genome sequence data
# Data

The DNA sequence data for this project are in `/uufs/chpc.utah.edu/common/home/gompert-group4/data/t_chumash_wgs/`. There are two batches of data: batch 1 = 192 samples, batch 2 = 288 samples (some of the latter are not part of the mapping and selection experiment). See [Nosil et al. 2020](https://www.proquest.com/docview/2473271380?pq-origsite=gscholar&fromopenview=true&sourcetype=Scholarly%20Journals) ([Nosil2020.pdf](https://github.com/user-attachments/files/17368867/Nosil2020.pdf)) for our analyses of GBS data from this same experiment.

**batch 1:** These samples were sequenced by GeT (Genotoul) and delivered July 2024. The IDs are at the start of the file names (Tchum-XXX, where XXX denotes the digits after 19_ in the actual file names). There are two files each of forward and reverse reads for each stick insect. 

**batch 2:** These samples were sequenced by GeT (Genotoul) and delivered October 2024. The IDs are again at the start of the file name but here are either 19_XXX (the actual ID we use) or GR806-XX for the non-experimental stick insects.

The experimental data (color and survival) are in the 2019_Tchumash* files in `/uufs/chpc.utah.edu/common/home/gompert-group4/data/t_chumash_wgs/`.

I am using our 2020 *T. chumash* genome assembly. This genome is described in [Nosil et al. 2020](https://www.proquest.com/docview/2473271380?pq-origsite=gscholar&fromopenview=true&sourcetype=Scholarly%20Journals). See [timema_chumash_29Feb2020_N4ago.report.pdf](https://github.com/user-attachments/files/17368853/timema_chumash_29Feb2020_N4ago.report.pdf) for a report from Dovetail. Past comparative alignments suggest general chromosome-level synteny with *T. crstinae* except that *T. cristinae* chromsomes 1, 3, 4 and 9 appear to be fused in *T. chumash*. Patrik suggested this is consistent with cytology (add source). Also note (from the report) that the genome size is notably larger than our other *Timema* genomes (2.4 Gbps) and the final assembly has a non-trivial number of duplicated BUSCOs (29). This could reflect some difficulities in the assembly (especially if the two genome copies differed enough) that might be cleaned up with phased genomes later. However, chromosome 8 appears to be largely syntenic (though not colinear) between *T. cristinae* and *T. chumash* and is the main focus of our interest ([see the TimemaFusion repository for past work on this](https://github.com/zgompert/TimemaFusion)).

# DNA sequence alignment

I aligned the whole genome sequence data to the *T. chumash* genome using `bwa-mem2` (version 2.0pre2).

I first indexed the genome.

```{bash}
/uufs/chpc.utah.edu/common/home/u6000989/source/bwa-mem2-2.0pre2_x64-linux/bwa-mem2 index mod_timema_chumash_29Feb2020_N4ago.fasta
```
I then ran the alignments. I copied all of the fastq files to scratch space for this.

```{bash}
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

This runs the following, which includes pipping the results on to `samtools` (version 1.16) to compress, sort and index the alignments.

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
        $fq1 =~ m/clean\/([A-Za-z0-9]+)/ or die "failed to match id $fq1\n";
        $ind = $1;
        $fq1 =~ m/([A-Za-z_\-0-9]+)_1\.fq\.gz$/ or die "failed match for file $fq1\n";
        $file = $1;
        system "/uufs/chpc.utah.edu/common/home/u6000989/source/bwa-mem2-2.0pre2_x64-linux/bwa-mem2 mem -t 1 -k 19 -r 1.5 -R \'\@RG\\tID:Tchum-"."$ind\\tLB:Tchum-"."$ind\\tSM:Tchum-"."$ind"."\' $genome $fq1 $fq2 | samtools sort -@ 2 -O BAM -o $ind"."_$file.bam - && samtools index -@ 2 $ind"."_$file.bam\n";

        $pm->finish;
}
```
