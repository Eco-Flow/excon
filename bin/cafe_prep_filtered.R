#!/usr/bin/Rscript

# CAFE preparation script with differential filtering
# This version is used when CAFE fails due to large size differentials

library(ape)
library(data.table)

# Get threshold and optional scale factor from command line
# args[1] = max_differential threshold (default 50)
# args[2] = branch-length scale factor applied after chronos() (default 1000)
args <- commandArgs(trailingOnly = TRUE)
max_differential <- if (length(args) >= 1) {
  as.numeric(args[1])
} else {
  as.numeric(Sys.getenv("CAFE_MAX_DIFF", "50"))
}
scale_factor <- if (length(args) >= 2) as.numeric(args[2]) else 1000

cat("================================================\n")
cat("CAFE PREP with Differential Filtering\n")
cat("Threshold:", max_differential, "\n")
cat("================================================\n\n")

tre <- read.tree('pruned_tree')
stopifnot(is.binary(tre))
stopifnot(is.rooted(tre))

if (!is.ultrametric(tre)) {
  tre <- chronoMPL(tre, SE = FALSE, test = FALSE)
}
tre$edge.length <- tre$edge.length * scale_factor
write.tree(tre, 'SpeciesTree_rooted_ultra.txt')

hog <- fread('N0.tsv')

# Handle both OrthoFinder v2 (HOG) and v3 (Orthogroup) column names
id_col <- if ('HOG' %in% names(hog)) 'HOG' else 'Orthogroup'

# Remove columns that may not exist in v3
if ('OG' %in% names(hog)) hog[, OG := NULL]
if ('Gene Tree Parent Clade' %in% names(hog)) hog[, `Gene Tree Parent Clade` := NULL]

hog <- melt(hog, id.vars=id_col, variable.name='species', value.name='pid')
hog <- hog[pid != '']
hog$n <- sapply(hog$pid, function(x) length(strsplit(x, ', ')[[1]]))

cat("Initial HOG count:", length(unique(hog[[id_col]])), "\n")

# Exclude HOGs with lots of genes in one or more species
keep <- hog[, list(n_max=max(n)), by=id_col][n_max < 100][[id_col]]
hog <- hog[get(id_col) %in% keep]
cat("After max<100 filter:", length(unique(hog[[id_col]])), "\n")

# Exclude HOGs present in only 1 species
keep <- hog[, .N, by=id_col][N > 1][[id_col]]
hog <- hog[get(id_col) %in% keep]
cat("After single-species filter:", length(unique(hog[[id_col]])), "\n")

# Calculate size differential (max - min) and filter
size_stats <- hog[, list(
  n_max = max(n),
  n_min = min(n),
  n_mean = mean(n),
  differential = max(n) - min(n)
), by = id_col]

# Report statistics before filtering
cat("\n--- Size Differential Statistics ---\n")
cat("Min differential:", min(size_stats$differential), "\n")
cat("Max differential:", max(size_stats$differential), "\n")
cat("Mean differential:", round(mean(size_stats$differential), 2), "\n")
cat("Median differential:", median(size_stats$differential), "\n")

# Show top problematic families
cat("\nTop 10 families with largest size differentials:\n")
top10 <- size_stats[order(-differential)][1:min(10, nrow(size_stats))]
print(top10[, .SD, .SDcols = c(id_col, 'n_min', 'n_max', 'differential')])

# Count families above threshold
n_above_threshold <- sum(size_stats$differential > max_differential)
cat("\nHOGs with differential >", max_differential, ":", n_above_threshold, "\n")

# Keep only HOGs with differential <= threshold
keep_diff <- size_stats[differential <= max_differential][[id_col]]
hog_filtered <- hog[get(id_col) %in% keep_diff]

cat("After differential filter:", length(unique(hog_filtered[[id_col]])), "\n")
cat("Families removed:", length(unique(hog[[id_col]])) - length(unique(hog_filtered[[id_col]])), "\n")
cat("Retention rate:", round(100 * length(unique(hog_filtered[[id_col]])) / length(unique(hog[[id_col]])), 2), "%\n\n")

# Write detailed filtering report
filtering_report <- size_stats[order(-differential)]
fwrite(filtering_report, 'hog_filtering_report.tsv', sep='\t')
cat("Detailed filtering report written to: hog_filtering_report.tsv\n")

# Create final count table — CAFE expects HOG as first column name
counts <- dcast(hog_filtered, get(id_col) ~ species, value.var='n', fill=0)
setnames(counts, 'id_col', 'HOG')
counts[, Desc := 'n/a']
setcolorder(counts, 'Desc')
fwrite(counts, 'hog_gene_counts.tsv', sep='\t')

cat("Filtered count table written to: hog_gene_counts.tsv\n")
cat("Final gene family count:", nrow(counts), "\n")

# Write the removed high-differential families for fixed-lambda re-analysis
large_ids <- size_stats[differential > max_differential][[id_col]]
if (length(large_ids) > 0) {
  hog_large <- hog[get(id_col) %in% large_ids]
  counts_large <- dcast(hog_large, get(id_col) ~ species, value.var='n', fill=0)
  setnames(counts_large, 'id_col', 'HOG')
  counts_large[, Desc := 'n/a']
  setcolorder(counts_large, 'Desc')
  fwrite(counts_large, 'hog_gene_counts_large.tsv', sep='\t')
  cat("Large-differential families written to: hog_gene_counts_large.tsv\n")
  cat("Count:", length(large_ids), "\n")
}

cat("================================================\n")
