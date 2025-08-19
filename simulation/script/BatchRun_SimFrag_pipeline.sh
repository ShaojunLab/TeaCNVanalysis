#!/bin/bash
#conda activate simfragenv


SCRIPT="./TeaCNV/simulation/script"
INPUT_DIR="./TeaCNV/simulation/simData"
OUTPUT_DIR="./TeaCNV/simulation/simData/Fragment"
[ -d $OUTPUT_DIR ] || mkdir -p $OUTPUT_DIR

clonal_groups=("Biclonal" "Triclonal" "Tetraclonal" "Pentaclonal")

tumor_sizes=(50)
for ((size=100; size<=1000; size+=100)); do
    tumor_sizes+=($size)
done
tumor_sizes+=(2000 5000 10000)

for group in "${clonal_groups[@]}"; do
  for size in "${tumor_sizes[@]}"; do
    for i in {1..10}; do
      current_path="${INPUT_DIR}/${group}/size${size}"
      out_path="${OUTPUT_DIR}/${group}/size${size}/${i}"
      [ -d ${out_path} ] || mkdir -p ${out_path}

      RDS1="${INPUT_DIR}/refCells/counts_ref${i}.rds"
      RDS2="${current_path}/${group}_counts_size${size}_${i}.rds"
      MERGED="${out_path}/merged_matrix.rds"
      OUTPUT="${out_path}/merged_fragments.tsv.gz"

      # Step 1: Merge matrices using R
      echo "Merging matrices..."
      Rscript "$SCRIPT/merge_rds_matrices.R" "$RDS1" "$RDS2" "$MERGED"
      if [[ $? -ne 0 ]]; then
        echo "Error during matrix merging."
        exit 1
      fi

      # Step 2: Simulate fragments using Python
      echo "Simulating fragments..."
      python3 "$SCRIPT/simulate_fragments_rds.py" --rds "$MERGED" --output "$OUTPUT"

      if [[ $? -ne 0 ]]; then
        echo "Error during fragment simulation."
        exit 1
      fi
    done
  done
done

echo "All done!"


