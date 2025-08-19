#!/usr/bin/env python3

import os
import re
import argparse
import gzip
import pysam
import pandas as pd
import numpy as np
from tqdm import tqdm
from collections import defaultdict
from rpy2.robjects import r, pandas2ri
from rpy2.robjects.conversion import localconverter
import rpy2.robjects.packages as rpackages



# Load required R packages
base = rpackages.importr('base')
utils = rpackages.importr('utils')
Matrix = rpackages.importr('Matrix')

def parse_peak_name(peak):
    if '-' in peak and peak.count('-') == 2:
        parts = peak.split('-')
        return parts[0], int(parts[1]), int(parts[2])
    elif ':' in peak and '-' in peak:
        # Format: chr1:9906-10568
        chrom, pos = peak.split(':')
        start, end = pos.split('-')
        return chrom, int(start), int(end)
    elif '_' in peak and peak.count('_') == 2:
        parts = peak.split('_')
        return parts[0], int(parts[1]), int(parts[2])
    else:
        raise ValueError(f"Unrecognized peak format: {peak}")

def simulate_fragments(rds_path, output_path):
    print(f"Reading RDS matrix from {rds_path}...")
    r_matrix = base.readRDS(rds_path)

    # Extract data from dgCMatrix
    peaks = list(r.rownames(r_matrix))
    cells = list(r.colnames(r_matrix))

    # Get sparse matrix slots
    matrix_data = r_matrix.do_slot("x")
    matrix_indices = r_matrix.do_slot("i")
    matrix_indptr = r_matrix.do_slot("p")

    print(f"Matrix dimensions: {len(peaks)} peaks x {len(cells)} cells")

    fragments = []

    for col_idx in tqdm(range(len(cells)), desc="Simulating fragments"):
        cell_barcode = cells[col_idx]
        start_ptr = matrix_indptr[col_idx]
        end_ptr = matrix_indptr[col_idx + 1]

        for data_idx in range(start_ptr, end_ptr):
            row_idx = matrix_indices[data_idx]
            count = int(matrix_data[data_idx])
            chrom, start, end = parse_peak_name(peaks[row_idx])

            for _ in range(count):
                insert_start = np.random.randint(start, end)
                insert_end = insert_start + np.random.randint(50, 500)  # Approx fragment length
                fragments.append((chrom, insert_start, insert_end, cell_barcode, 1))

    
    # 
    frag_df = pd.DataFrame(fragments, columns=["chr", "start", "end", "barcode", "score"])
    frag_df.sort_values(by=["chr", "start", "end"], inplace=True)

    print(f"Writing BGZF compressed fragments to {output_path} with header...")

    # pysam.BGZFile BGZF
    with pysam.BGZFile(output_path, "w") as fout:
        # header = "\t".join(["chr", "start", "end", "barcode", "score"]) + "\n"
        # fout.write(header.encode())
        for row in frag_df.itertuples(index=False):
            line = "\t".join(map(str, row)) + "\n"
            fout.write(line.encode())

    print("Indexing fragment file...")
    pysam.tabix_index(output_path, preset="bed", force=True)
    print("Done.")

    # Write summary statistics
    summary_path = output_path.replace(".tsv.gz", ".summary.tsv")
    print(f"Writing summary to {summary_path}...")
    frag_df = pd.DataFrame(fragments, columns=["chr", "start", "end", "barcode", "score"])
    stats = frag_df.groupby("barcode").size().reset_index(name="fragment_count")
    stats.to_csv(summary_path, sep="\t", index=False)
    print("Summary written.")


if __name__ == "__main__":
    parser = argparse.ArgumentParser(description="Simulate fragment file from sparse peak-cell RDS matrix")
    parser.add_argument("--rds", required=True, help="Path to .rds matrix file")
    parser.add_argument("--output", default="simulated_fragments.tsv.gz", help="Output fragment file name")
    args = parser.parse_args()

    try:
        simulate_fragments(args.rds, args.output)
    except Exception as e:
        print(f"Error during fragment simulation: {e}")

