#!/bin/bash 
#SBATCH --time=168:00:00
#SBATCH --nodes=1
#SBATCH --ntasks=96
#SBATCH --mem=600000
#SBATCH --account=gompert
#SBATCH --qos=gompert-grn
#SBATCH --partition=gompert-grn
#SBATCH --job-name=cactus
#SBATCH --mail-type=FAIL
#SBATCH --mail-user=zach.gompert@usu.edu

module load cactus

cd /scratch/general/nfs1/u6000989/cactus

# Max number of concurrent jobs
MAX_JOBS=4

# Input files and associated genome IDs
comps=(
"cactusTchum_E240159H1_E240161H1.txt"
"cactusTchum_E240159H1_E240161H2.txt"
"cactusTchum_E240159H1_E240163H1.txt"
"cactusTchum_E240159H1_E240163H2.txt"
"cactusTchum_E240159H2_E240161H1.txt"
"cactusTchum_E240159H2_E240161H2.txt"
"cactusTchum_E240159H2_E240163H1.txt"
"cactusTchum_E240159H2_E240163H2.txt"
"cactusTchum_E240161H1_E240163H1.txt"
"cactusTchum_E240161H1_E240163H2.txt"
"cactusTchum_E240161H2_E240163H1.txt"
"cactusTchum_E240161H2_E240163H2.txt"
)

ids1=(
  "t_chum_24_0159_h1"
  "t_chum_24_0159_h1"
  "t_chum_24_0159_h1"
  "t_chum_24_0159_h1"
  "t_chum_24_0159_h2"
  "t_chum_24_0159_h2"
  "t_chum_24_0159_h2"
  "t_chum_24_0159_h2"
  "t_chum_24_0161_h1"
  "t_chum_24_0161_h1"
  "t_chum_24_0161_h2"
  "t_chum_24_0161_h2"
)

ids2=(
  "t_chum_24_0161_h1"
  "t_chum_24_0161_h2"
  "t_chum_24_0163_h1"
  "t_chum_24_0163_h2"
  "t_chum_24_0161_h1"
  "t_chum_24_0161_h2"
  "t_chum_24_0163_h1"
  "t_chum_24_0163_h2"
  "t_chum_24_0163_h1"
  "t_chum_24_0163_h2"
  "t_chum_24_0163_h1"
  "t_chum_24_0163_h2"
)

# Loop over jobs
for i in "${!comps[@]}"; do
  file="${comps[$i]}"
  id1="${ids1[$i]}"
  id2="${ids2[$i]}"

  if [[ "$file" =~ cactusTchum_([a-zA-Z0-9]+)_([a-zA-Z0-9_]+)\.txt ]]; then
    base="${BASH_REMATCH[1]}_${BASH_REMATCH[2]}"
  else
    echo "Failed to parse: $file"
    continue
  fi

  (
    cactus "jobStore_$base" \
      "/uufs/chpc.utah.edu/common/home/gompert-group4/data/timema/hic_genomes/comp_aligns/$file" \
      "cactus${base}.hal" --maxCores 24

    ~/source/hal/bin/halSynteny --queryGenome "$id2" --targetGenome "$id1" \
      "cactus${base}.hal" "out_synteny_${base}.psl"
  ) &

  # Limit the number of background jobs
  while (( $(jobs -rp | wc -l) >= MAX_JOBS )); do
    wait -n
  done
done

# Wait for all remaining background jobs to finish
wait


