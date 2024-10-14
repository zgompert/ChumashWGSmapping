# ChumashWGSmapping
Genomic analyses of color and survival in a selection experiment for Timema chumash based on whole genome sequence data
# Data

The DNA sequence data for this project are in `/uufs/chpc.utah.edu/common/home/gompert-group4/data/t_chumash_wgs/`. There are two batches of data: batch 1 = 192 samples, batch 2 = 288 samples (some of the latter are not part of the mapping and selection experiment). See [Nosil et al. 2020](https://www.proquest.com/docview/2473271380?pq-origsite=gscholar&fromopenview=true&sourcetype=Scholarly%20Journals) ([Nosil2020.pdf](https://github.com/user-attachments/files/17368867/Nosil2020.pdf)) for our analyses of GBS data from this same experiment.

**batch 1:** These samples were sequenced by GeT (Genotoul) and delivered July 2024. The IDs are at the start of the file names (Tchum-XXX, where XXX denotes the digits after 19_ in the actual file names). There are two files each of forward and reverse reads for each stick insect. 

**batch 2:** These samples were sequenced by GeT (Genotoul) and delivered October 2024. The IDs are again at the start of the file name but here are either 19_XXX (the actual ID we use) or GR806-XX for the non-experimental stick insects.

The experimental data (color and survival) are in the 2019_Tchumash* files in `/uufs/chpc.utah.edu/common/home/gompert-group4/data/t_chumash_wgs/`.

I am using our 2020 *T. chumash* genome assembly. This genome is described in [Nosil et al. 2020](https://www.proquest.com/docview/2473271380?pq-origsite=gscholar&fromopenview=true&sourcetype=Scholarly%20Journals). See [timema_chumash_29Feb2020_N4ago.report.pdf](https://github.com/user-attachments/files/17368853/timema_chumash_29Feb2020_N4ago.report.pdf) for a report from Dovetail. Past comparative alignments suggest general chromosome-level synteny with *T. crstinae* except that *T. cristinae* chromsomes 1, 3, 4 and 9 appear to be fused in *T. chumash*. Patrik suggested this is consistent with cytology (add source). Also note (from the report) that the genome size is notably larger than our other *Timema* genomes (2.4 Gbps) and the final assembly has a non-trivial number of duplicated BUSCOs (29). This could reflect some difficulities in the assembly (especially if the two genome copies differed enough) that might be cleaned up with phased genomes later. However, chromosome 8 appears to be largely syntenic (though not colinear) between *T. cristinae* and *T. chumash* and is the main focus of our interest ([see the TimemaFusion repository for past work on this](https://github.com/zgompert/TimemaFusion)).

# DNA sequence alignment

I aligned the whole genome sequence data to the *T. chumash* genome using `bwa-mem2`.


