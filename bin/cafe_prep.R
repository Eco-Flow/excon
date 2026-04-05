#!/usr/bin/Rscript
library(ape)
library(data.table)

# Optional first arg: branch-length scale factor (default 1000).
# Applied *after* chronos() so the scaling is not lost when the tree is
# non-ultrametric (chronos recalculates branch lengths from scratch).
args <- commandArgs(trailingOnly = TRUE)
scale_factor <- if (length(args) >= 1) as.numeric(args[1]) else 1000

tre <- read.tree('pruned_tree')
stopifnot(is.binary(tre))
stopifnot(is.rooted(tre))
if (!is.ultrametric(tre)) {
  tre <- chronos(tre)
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

# Exclude HOGs with lots of genes in one or more species
keep <- hog[, list(n_max=max(n)), by=id_col][n_max < 100][[id_col]]
hog <- hog[get(id_col) %in% keep]

# Exclude HOGs present in only 1 species
keep <- hog[, .N, by=id_col][N > 1][[id_col]]
hog <- hog[get(id_col) %in% keep]

counts <- dcast(hog, get(id_col) ~ species, value.var='n', fill=0)
setnames(counts, 'id_col', 'HOG')  # CAFE expects HOG as the first column name
counts[, Desc := 'n/a']
setcolorder(counts, 'Desc')
fwrite(counts, 'hog_gene_counts.tsv', sep='\t')
