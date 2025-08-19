#!/bin/bash
#conda activate simfragenv


SCRIPT="./TeaCNV/simulation/script"
INPUT_DIR="./TeaCNV/simulation/simData_diffCNVpct"
OUTPUT_DIR="./TeaCNV/simulation/simData_diffCNVpct/Fragment"
[ -d $OUTPUT_DIR ] || mkdir -p $OUTPUT_DIR

clonal_groups=("Monoclonal" "Biclonal" "Triclonal" "Tetraclonal")
cnvPct=(0.3 0.4 0.6)

for pct in "${cnvPct[@]}"; do
  for group in "${clonal_groups[@]}"; do
    current_path="${INPUT_DIR}/CNVpct${pct}/${group}"
    for i in {1..10}; do  
      out_path="${OUTPUT_DIR}/CNVpct${pct}/${group}/${i}"
      [ -d ${out_path} ] || mkdir -p ${out_path}

      RDS1="${INPUT_DIR}/refCells/counts_ref${i}.rds"
      RDS2="${current_path}/counts_${i}.rds"
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

