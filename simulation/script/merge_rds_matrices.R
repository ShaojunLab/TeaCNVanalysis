#!/usr/bin/env Rscript

args <- commandArgs(trailingOnly = TRUE)
if (length(args) != 3) {
  stop("Usage: Rscript merge_rds_matrices.R <input1.rds> <input2.rds> <output.rds>")
}

input1 <- args[1]
input2 <- args[2]
output <- args[3]
library(Matrix)

cat("Reading input matrices...\n")
mat1 <- readRDS(input1)
mat2 <- readRDS(input2)

# Check rownames match (peaks)
if (!identical(as.character(rownames(mat1)), as.character(rownames(mat2)))) {
  stop("Peak rownames do not match.")
}

# Add prefixes to barcodes to keep unique identity
colnames(mat1) <- paste0("sample1_", colnames(mat1))
colnames(mat2) <- paste0("sample2_", colnames(mat2))

cat("Merging matrices by columns...\n")
merged <- cbind(mat1, mat2)

cat("Saving merged matrix to", output, "...\n")
saveRDS(merged, file = output)

cat("Merge complete.\n")
