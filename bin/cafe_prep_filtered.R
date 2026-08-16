#!/usr/bin/Rscript

# CAFE input preparation WITH differential filtering.
# Used on retry when the un-filtered base run fails on a large size differential.
# Mirrors cafe_prep.R (subsetting, dated-tree handling, reasoned report) and adds
# a max-minus-min differential filter, routing over-threshold families to a
# fixed-lambda re-analysis (hog_gene_counts_large.tsv).
#
# Positional args (supplied by cafe_prep.nf):
#   args[1] = max_differential threshold (max-min copies) (default 50)
#   args[2] = branch-length scale factor applied AFTER chronoMPL() (default 1000)
#             Ignored when the tree is already dated (args[3] == "true").
#   args[3] = "true" if --input_tree_is_dated (tree validated + used unchanged)

library(ape)
library(data.table)

args <- commandArgs(trailingOnly = TRUE)
max_differential <- if (length(args) >= 1) as.numeric(args[1]) else
  as.numeric(Sys.getenv("CAFE_MAX_DIFF", "50"))
scale_factor <- if (length(args) >= 2) as.numeric(args[2]) else 1000
is_dated     <- if (length(args) >= 3) tolower(args[3]) %in% c("true", "1", "yes") else FALSE

cat("================================================\n")
cat("CAFE PREP with Differential Filtering\n")
cat("Threshold (max-min):", max_differential, "\n")
cat("Dated input tree:", is_dated, "\n")
cat("================================================\n\n")

## ------------------------------------------------------------------
## 1. Tree
## ------------------------------------------------------------------
tre <- read.tree('pruned_tree')
stopifnot("Input tree must be rooted" = is.rooted(tre))
stopifnot("Input tree must be binary" = is.binary(tre))
if (any(is.na(tre$edge.length)) || any(tre$edge.length <= 0)) {
  stop("Input tree has missing or non-positive branch lengths.")
}

if (is_dated) {
  if (!is.ultrametric(tre, tol = 1e-3)) {
    stop("--input_tree_is_dated set but supplied tree is not ultrametric.")
  }
} else {
  if (!is.ultrametric(tre)) tre <- chronoMPL(tre)
  tre$edge.length <- tre$edge.length * scale_factor
}
tree_leaves <- tre$tip.label

## ------------------------------------------------------------------
## 2. N0 subset to tree leaves
## ------------------------------------------------------------------
hog_wide <- fread('N0.tsv')
id_col <- if ('HOG' %in% names(hog_wide)) 'HOG' else 'Orthogroup'
for (drop_col in c('OG', 'Gene Tree Parent Clade')) {
  if (drop_col %in% names(hog_wide)) hog_wide[, (drop_col) := NULL]
}
species_cols <- setdiff(names(hog_wide), id_col)
missing_leaves <- setdiff(tree_leaves, species_cols)
if (length(missing_leaves) > 0) {
  stop("Tree leaves absent from N0.tsv columns: ",
       paste(missing_leaves, collapse = ", "))
}
hog_wide <- hog_wide[, c(id_col, tree_leaves), with = FALSE]
all_hogs <- hog_wide[[id_col]]

## ------------------------------------------------------------------
## 3. Long form -> complete count matrix
## ------------------------------------------------------------------
hog <- melt(hog_wide, id.vars = id_col, variable.name = 'species', value.name = 'pid')
hog <- hog[pid != '']
hog[, n := vapply(strsplit(pid, ', '), length, integer(1))]
hog[, species := factor(species, levels = tree_leaves)]
counts <- dcast(hog, get(id_col) ~ species, value.var = 'n', fill = 0, drop = c(TRUE, FALSE))
setnames(counts, 'id_col', 'HOG')
for (leaf in tree_leaves) if (!leaf %in% names(counts)) counts[, (leaf) := 0L]
setcolorder(counts, c('HOG', tree_leaves))

## ------------------------------------------------------------------
## 4. Stats + exclusion reasons (incl. differential threshold)
## ------------------------------------------------------------------
cmat <- as.matrix(counts[, ..tree_leaves])
stats <- data.table(
  HOG               = counts$HOG,
  n_species_present = rowSums(cmat > 0),
  n_max             = apply(cmat, 1, max),
  n_min             = apply(cmat, 1, min),
  total_genes       = rowSums(cmat)
)
stats[, differential := n_max - n_min]

empty_hogs <- setdiff(all_hogs, counts$HOG)
if (length(empty_hogs) > 0) {
  stats <- rbind(stats, data.table(
    HOG = empty_hogs, n_species_present = 0L, n_max = 0L, n_min = 0L,
    total_genes = 0L, differential = 0L
  ))
}

MAX_COPIES <- 100L
reason <- rep('retained', nrow(stats))
reason[stats$differential > max_differential]  <- 'differential_gt_threshold'
reason[stats$n_max >= MAX_COPIES]              <- 'max_copies_ge_100'
reason[stats$n_species_present == 1]           <- 'single_species'
reason[stats$n_species_present == 0]           <- 'empty_after_subset'
stats[, exclusion_reason := reason]
stats[, excluded := exclusion_reason != 'retained']

setorder(stats, exclusion_reason, -differential)
fwrite(stats, 'hog_filtering_report.tsv', sep = '\t')
cat("\n--- Filtering summary (reason : n HOGs) ---\n")
print(stats[, .N, by = exclusion_reason])

## ------------------------------------------------------------------
## 5. Retained counts + tree + large-family table
## ------------------------------------------------------------------
keep_hogs <- stats[exclusion_reason == 'retained', HOG]
counts_keep <- counts[HOG %in% keep_hogs]
counts_keep[, Desc := 'n/a']
setcolorder(counts_keep, c('Desc', 'HOG', tree_leaves))
fwrite(counts_keep, 'hog_gene_counts.tsv', sep = '\t')
cat("Retained families:", nrow(counts_keep), "\n")

stopifnot("Tree tips must match gene-count columns" =
            setequal(tre$tip.label, tree_leaves))
write.tree(tre, 'cafe_input_tree.txt')
write.tree(tre, 'SpeciesTree_rooted_ultra.txt')

# High-differential families for fixed-lambda re-analysis
large_ids <- stats[exclusion_reason == 'differential_gt_threshold', HOG]
if (length(large_ids) > 0) {
  counts_large <- counts[HOG %in% large_ids]
  counts_large[, Desc := 'n/a']
  setcolorder(counts_large, c('Desc', 'HOG', tree_leaves))
  fwrite(counts_large, 'hog_gene_counts_large.tsv', sep = '\t')
  cat("Large-differential families written to hog_gene_counts_large.tsv:",
      length(large_ids), "\n")
}
cat("================================================\n")
