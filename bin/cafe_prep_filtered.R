#!/usr/bin/Rscript

# CAFE preparation script with differential filtering
# This version is used when CAFE fails due to large size differentials

library(ape)
library(data.table)

# Get threshold from command line or environment variable
args <- commandArgs(trailingOnly = TRUE)
max_differential <- if (length(args) >= 1) {
    as.numeric(args[1])
} else {
    as.numeric(Sys.getenv("CAFE_MAX_DIFF", "50"))
}

cat("================================================\n")
cat("CAFE PREP with Differential Filtering\n")
cat("Threshold:", max_differential, "\n")
cat("================================================\n\n")

tre <- read.tree('pruned_tree')
stopifnot(is.binary(tre))
stopifnot(is.rooted(tre))

if(is.ultrametric(tre)) {
    utre <- tre
} else{
    utre <- chronos(tre)
}
write.tree(utre, 'SpeciesTree_rooted_ultra.txt')

hog <- fread('N0.tsv')
hog[, OG := NULL]
hog[, `Gene Tree Parent Clade` := NULL]
hog <- melt(hog, id.vars='HOG', variable.name='species', value.name='pid')
hog <- hog[pid != '']
hog$n <- sapply(hog$pid, function(x) length(strsplit(x, ', ')[[1]]))

cat("Initial HOG count:", length(unique(hog$HOG)), "\n")

# Exclude HOGs with lots of genes in a one or more species. 
# See also cafe tutorial about filtering gene families
keep <- hog[, list(n_max=max(n)), HOG][n_max < 100]$HOG
hog <- hog[HOG %in% keep]
cat("After max<100 filter:", length(unique(hog$HOG)), "\n")

# Exclude HOGs present in only 1 species
keep <- hog[, .N, HOG][N > 1]$HOG
hog <- hog[HOG %in% keep]
cat("After single-species filter:", length(unique(hog$HOG)), "\n")

# NEW: Calculate size differential (max - min) and filter
# This addresses the CAFE5 warning about large size differentials
size_stats <- hog[, list(
    n_max = max(n),
    n_min = min(n),
    n_mean = mean(n),
    differential = max(n) - min(n)
), by = HOG]

# Report statistics before filtering
cat("\n--- Size Differential Statistics ---\n")
cat("Min differential:", min(size_stats$differential), "\n")
cat("Max differential:", max(size_stats$differential), "\n")
cat("Mean differential:", round(mean(size_stats$differential), 2), "\n")
cat("Median differential:", median(size_stats$differential), "\n")

# Show top problematic families
cat("\nTop 10 families with largest size differentials:\n")
top10 <- size_stats[order(-differential)][1:min(10, nrow(size_stats))]
print(top10[, .(HOG, n_min, n_max, differential)])

# Count families above threshold
n_above_threshold <- sum(size_stats$differential > max_differential)
cat("\nHOGs with differential >", max_differential, ":", n_above_threshold, "\n")

# Keep only HOGs with differential <= threshold
keep_diff <- size_stats[differential <= max_differential]$HOG
hog_filtered <- hog[HOG %in% keep_diff]

cat("After differential filter:", length(unique(hog_filtered$HOG)), "\n")
cat("Families removed:", length(unique(hog$HOG)) - length(unique(hog_filtered$HOG)), "\n")
cat("Retention rate:", 
    round(100 * length(unique(hog_filtered$HOG)) / length(unique(hog$HOG)), 2), 
    "%\n\n")

# Write detailed filtering report
filtering_report <- size_stats[order(-differential)]
fwrite(filtering_report, 'hog_filtering_report.tsv', sep='\t')
cat("✅ Detailed filtering report written to: hog_filtering_report.tsv\n")

# Create final count table
counts <- dcast(hog_filtered, HOG ~ species, value.var='n', fill=0)
counts[, Desc := 'n/a']
setcolorder(counts, 'Desc')
fwrite(counts, 'hog_gene_counts.tsv', sep='\t')

cat("✅ Filtered count table written to: hog_gene_counts.tsv\n")
cat("Final gene family count:", nrow(counts), "\n")
cat("================================================\n")
