#!/usr/bin/env Rscript
# cafe_large_families.R
# Gene-count table for the families CAFE_PREP left out of the main CAFE5 run for
# their size, so CAFE_RUN_LARGE can fit each one on its own:
#   max_copies_ge_100          (>= 100 copies in some species; any CAFE_PREP attempt)
#   differential_gt_threshold  (max - min copies over --cafe_max_differential; retries)
# Counts are built from N0.tsv exactly as cafe_prep.R does.
#
# Usage: cafe_large_families.R <N0.tsv> <hog_filtering_report.tsv> <cafe_input_tree.txt>
# Writes hog_gene_counts_large.tsv when there is at least one such family, and
# large_families.txt (one line, for the run log) either way.

suppressPackageStartupMessages({
  library(ape)
  library(data.table)
})

args   <- commandArgs(trailingOnly = TRUE)
leaves <- read.tree(args[3])$tip.label
report <- fread(args[2])
sizes  <- c("max_copies_ge_100", "differential_gt_threshold")
large  <- report[exclusion_reason %in% sizes]

if (nrow(large) == 0) {
  writeLines("Large-family track: no family was left out of the main CAFE5 run for its size, so there is nothing to run.",
             "large_families.txt")
  quit(status = 0)
}

n0     <- fread(args[1])
id_col <- if ("HOG" %in% names(n0)) "HOG" else "Orthogroup"
n0     <- n0[get(id_col) %in% large$HOG, c(id_col, leaves), with = FALSE]

# One count per species: the number of genes listed in its cell
count_genes <- function(cell) ifelse(cell == "", 0L, lengths(strsplit(cell, ", ")))
counts <- n0[, c(list(HOG = get(id_col)), lapply(.SD, count_genes)), .SDcols = leaves]
counts[, Desc := "n/a"]
setcolorder(counts, c("Desc", "HOG", leaves))
fwrite(counts, "hog_gene_counts_large.tsv", sep = "\t")

by_reason <- large[, .N, by = exclusion_reason]
writeLines(sprintf("Large-family track: %d families left out of the main CAFE5 run for their size (%s) are fitted one at a time.",
                   nrow(counts), paste(sprintf("%d %s", by_reason$N, by_reason$exclusion_reason), collapse = ", ")),
           "large_families.txt")
