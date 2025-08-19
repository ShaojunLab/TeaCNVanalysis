#!/usr/bin/env python3

import argparse
import gzip
import random
from scipy.io import mmread
from scipy import sparse

def read_peaks(bed_file):
    peaks = []
    open_fn = gzip.open if bed_file.endswith(".gz") else open
    with open_fn(bed_file, "rt") as f:
        for line in f:
            chrom, start, end = line.strip().split()[:3]
            peaks.append((chrom, int(start), int(end)))
    return peaks

def read_barcodes(barcode_file):
    barcodes = []
    open_fn = gzip.open if barcode_file.endswith(".gz") else open
    with open_fn(barcode_file, "rt") as f:
        barcodes = [line.strip() for line in f]
    return barcodes

def simulate_fragments(peaks, barcodes, matrix, output_file):
    with gzip.open(output_file, "wt") as out:
        for i, j, v in zip(matrix.row, matrix.col, matrix.data):
            chrom, start, end = peaks[i]
            barcode = barcodes[j]
            for _ in range(int(v)):
                frag_start = random.randint(start, end - 25) if (end - start > 25) else start
                frag_len = random.randint(25, 150)
                frag_end = min(frag_start + frag_len, end)
                strand = random.choice(["+", "-"])
                out.write(f"{chrom}\t{frag_start}\t{frag_end}\t{barcode}\t1\t{strand}\n")

def main():
    parser = argparse.ArgumentParser(description="Simulate a fragments.tsv.gz file from CellRanger ATAC count matrix")
    parser.add_argument("--matrix", required=True, help="Path to matrix.mtx.gz")
    parser.add_argument("--barcodes", required=True, help="Path to barcodes.tsv.gz")
    parser.add_argument("--peaks", required=True, help="Path to peaks.bed (or .gz)")
    parser.add_argument("--output", default="simulated_fragments.tsv.gz", help="Output fragments file (default: simulated_fragments.tsv.gz)")
    args = parser.parse_args()

    print("Loading matrix...")
    matrix = mmread(args.matrix).tocoo()

    print("Loading barcodes...")
    barcodes = read_barcodes(args.barcodes)

    print("Loading peaks...")
    peaks = read_peaks(args.peaks)

    assert matrix.shape[0] == len(peaks), "Matrix row count does not match number of peaks"
    assert matrix.shape[1] == len(barcodes), "Matrix column count does not match number of barcodes"

    print("Generating simulated fragments...")
    simulate_fragments(peaks, barcodes, matrix, args.output)

    print(f"Done! Simulated fragments saved to: {args.output}")

if __name__ == "__main__":
    main()
