#!/bin/bash
#conda activate simfragenv


SCRIPT="./TeaCNV/simulation/script"
INPUT_DIR="./TeaCNV/simulation/simData"
OUTPUT_DIR="./TeaCNV/simulation/simData/Fragment"
[ -d $OUTPUT_DIR ] || mkdir -p $OUTPUT_DIR

RDS1=$INPUT_DIR/refCells/counts_ref1.rds
RDS2=$INPUT_DIR/Biclonal/size1000/Biclonal_counts_size1000_1.rds
MERGED=$OUTPUT_DIR/merged_matrix.rds
OUTPUT=$OUTPUT_DIR/merged_fragments.tsv.gz

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

echo "All done!"


