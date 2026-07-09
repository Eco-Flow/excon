#!/usr/bin/Rscript

# CAFE input preparation.
#
# Responsibilities:
#   1. Load the species tree that CAFE will use and validate it.
#   2. Subset the OrthoFinder N0.tsv gene-count columns to exactly the species
#      present in that tree, retaining the original HOG identifiers.
#   3. Filter families (empty-after-subset -> single-species -> very large),
#      recording the reason for every exclusion in hog_filtering_report.tsv.
#   4. Emit the gene-count table (hog_gene_counts.tsv) and the exact tree CAFE
#      should consume (cafe_input_tree.txt), whose leaf set matches the counts.
#
# Positional args (all optional, supplied by cafe_prep.nf):
#   args[1] = branch-length scale factor applied AFTER chronoMPL() (default 1000).
#             Ignored entirely when the tree is already dated (args[2] == "true").
#   args[2] = "true" if --input_tree_is_dated: the tree is already time-calibrated
#             and ultrametric; it is validated and written out UNCHANGED
#             (no chronoMPL(), no rescaling, effectively scale factor 1).

library(ape)
library(data.table)

args <- commandArgs(trailingOnly = TRUE)
scale_factor <- if (length(args) >= 1) as.numeric(args[1]) else 1000
is_dated     <- if (length(args) >= 2) tolower(args[2]) %in% c("true", "1", "yes") else FALSE

cat("================================================\n")
cat("CAFE PREP\n")
cat("Dated (ultrametric, time-calibrated) input tree:", is_dated, "\n")
if (!is_dated) cat("Post-chronoMPL scale factor:", scale_factor, "\n")
cat("================================================\n\n")

## ------------------------------------------------------------------
## 1. Tree: validate, and produce the exact tree CAFE will consume.
## ------------------------------------------------------------------
tre <- read.tree('pruned_tree')

stopifnot("Input tree must be rooted"  = is.rooted(tre))
stopifnot("Input tree must be binary"  = is.binary(tre))
if (any(is.na(tre$edge.length)) || any(tre$edge.length <= 0)) {
  stop("Input tree has missing or non-positive branch lengths — CAFE requires all branches > 0.")
}

if (is_dated) {
  # Trust the supplied calibration: validate ultrametricity, change nothing.
  tip_depths <- node.depth.edgelength(tre)[1:Ntip(tre)]
  if (!is.ultrametric(tre, tol = 1e-3)) {
    stop("--input_tree_is_dated was set but the supplied tree is NOT ultrametric ",
         "(root-to-tip depth range ", signif(min(tip_depths), 4), " to ",
         signif(max(tip_depths), 4),
         "). Supply a properly time-calibrated ultrametric tree or drop the flag.")
  }
  cat("Dated tree accepted unchanged (no chronoMPL, no rescaling).\n")
} else {
  # Legacy behaviour: make ultrametric via mean-path-length, then scale ONCE.
  if (!is.ultrametric(tre)) {
    cat("Tree is not ultrametric — applying chronoMPL().\n")
    tre <- chronoMPL(tre)
  }
  tre$edge.length <- tre$edge.length * scale_factor
}

tree_leaves <- tre$tip.label
cat("Tree leaves (", length(tree_leaves), "):\n", paste(tree_leaves, collapse = ", "), "\n\n", sep = "")

## ------------------------------------------------------------------
## 2. Load N0 and subset species columns to the tree leaves.
## ------------------------------------------------------------------
hog_wide <- fread('N0.tsv')

# Handle both OrthoFinder v2 (HOG) and v3 (Orthogroup) column names
id_col <- if ('HOG' %in% names(hog_wide)) 'HOG' else 'Orthogroup'

# Remove non-species metadata columns that may or may not be present
for (drop_col in c('OG', 'Gene Tree Parent Clade')) {
  if (drop_col %in% names(hog_wide)) hog_wide[, (drop_col) := NULL]
}

species_cols <- setdiff(names(hog_wide), id_col)

# Every tree leaf MUST have a matching gene-count column.
missing_leaves <- setdiff(tree_leaves, species_cols)
if (length(missing_leaves) > 0) {
  stop("The following tree leaves have no matching column in N0.tsv:\n  ",
       paste(missing_leaves, collapse = ", "),
       "\nN0.tsv species columns are:\n  ", paste(species_cols, collapse = ", "),
       "\nLeaf names and gene-count column headers must match exactly.")
}

dropped_species <- setdiff(species_cols, tree_leaves)
if (length(dropped_species) > 0) {
  cat("Subsetting N0.tsv to the", length(tree_leaves), "tree species; dropping",
      length(dropped_species), "column(s) not in the tree:\n  ",
      paste(dropped_species, collapse = ", "), "\n\n")
}

# Keep only the id column + the tree's species columns
hog_wide <- hog_wide[, c(id_col, tree_leaves), with = FALSE]

all_hogs <- hog_wide[[id_col]]
cat("HOGs in (subset) N0.tsv:", length(all_hogs), "\n")

## ------------------------------------------------------------------
## 3. Long form -> per-leaf gene counts -> complete count matrix.
## ------------------------------------------------------------------
hog <- melt(hog_wide, id.vars = id_col, variable.name = 'species', value.name = 'pid')
hog <- hog[pid != '']
hog[, n := vapply(strsplit(pid, ', '), length, integer(1))]

# Wide count matrix (rows = non-empty HOGs, cols = every tree leaf, fill 0)
hog[, species := factor(species, levels = tree_leaves)]
counts <- dcast(hog, get(id_col) ~ species, value.var = 'n', fill = 0, drop = c(TRUE, FALSE))
setnames(counts, 'id_col', 'HOG')

# Guarantee a column for every tree leaf (a leaf absent from all HOGs gets a 0 column)
for (leaf in tree_leaves) if (!leaf %in% names(counts)) counts[, (leaf) := 0L]
setcolorder(counts, c('HOG', tree_leaves))

## ------------------------------------------------------------------
## 4. Per-HOG statistics and exclusion reasons.
## ------------------------------------------------------------------
cmat <- as.matrix(counts[, ..tree_leaves])
stats <- data.table(
  HOG               = counts$HOG,
  n_species_present = rowSums(cmat > 0),
  n_max             = apply(cmat, 1, max),
  n_min             = apply(cmat, 1, min),   # min over ALL leaves (0 if absent anywhere)
  total_genes       = rowSums(cmat)
)
stats[, differential := n_max - n_min]

# HOGs present in N0 but empty after subsetting never made it into `counts`.
empty_hogs <- setdiff(all_hogs, counts$HOG)
if (length(empty_hogs) > 0) {
  stats <- rbind(stats, data.table(
    HOG = empty_hogs, n_species_present = 0L, n_max = 0L, n_min = 0L,
    total_genes = 0L, differential = 0L
  ))
}

MAX_COPIES <- 100L   # families with >=100 copies in any one species are excluded

reason <- rep('retained', nrow(stats))
reason[stats$n_max >= MAX_COPIES]     <- 'max_copies_ge_100'
reason[stats$n_species_present == 1]  <- 'single_species'
reason[stats$n_species_present == 0]  <- 'empty_after_subset'
stats[, exclusion_reason := reason]
stats[, excluded := exclusion_reason != 'retained']

setorder(stats, exclusion_reason, -differential)
fwrite(stats, 'hog_filtering_report.tsv', sep = '\t')

cat("\n--- Filtering summary (reason : n HOGs) ---\n")
print(stats[, .N, by = exclusion_reason])
cat("Detailed per-HOG report written to: hog_filtering_report.tsv\n\n")

## ------------------------------------------------------------------
## 5. Write the retained count table and the CAFE tree.
## ------------------------------------------------------------------
keep_hogs <- stats[exclusion_reason == 'retained', HOG]
counts <- counts[HOG %in% keep_hogs]

counts[, Desc := 'n/a']
setcolorder(counts, c('Desc', 'HOG', tree_leaves))
fwrite(counts, 'hog_gene_counts.tsv', sep = '\t')
cat("Retained families:", nrow(counts), "\n")
cat("Gene-count table written to: hog_gene_counts.tsv\n")

# The tree CAFE actually consumes. Its tips match the count columns exactly.
stopifnot("Tree tips must match gene-count columns" =
            setequal(tre$tip.label, tree_leaves))
write.tree(tre, 'cafe_input_tree.txt')
# Backwards-compatible alias retained for downstream/publish consumers.
write.tree(tre, 'SpeciesTree_rooted_ultra.txt')
cat("CAFE input tree written to: cafe_input_tree.txt (also SpeciesTree_rooted_ultra.txt)\n")
cat("================================================\n")
