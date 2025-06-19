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

