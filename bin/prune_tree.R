#!/usr/bin/env Rscript

# Prune a species tree to a subset of taxa, so CAFE5 can be run on one clade at
# a time.
#
# Running CAFE on a subset is often preferable to running it across a whole
# analysis: distantly related outgroups leave many families empty at the root,
# and CAFE5 discards those, which tends to remove exactly the lineage-specific
# families of interest. cafe_prep.R subsets the gene-count table to whatever
# species the tree contains, so pruning the tree is all that is required.
#
# Pruning preserves ultrametricity and node ages, so a time-calibrated tree stays
# calibrated and the ages of retained nodes do not change.
#
# Arguments (positional):
#   args[1] = input tree (rooted Newick)
#   args[2] = selection mode: "clade" or "species"
#   args[3] = selection:
#             clade   - two or more tip names; everything descended from their
#                       most recent common ancestor is kept
#             species - a comma-separated list of tip names, or a file with one
#                       name per line; exactly those tips are kept
#
# Outputs: SpeciesTree_pruned.nwk, pruned_tree_species.tsv

suppressPackageStartupMessages(library(ape))

args      <- commandArgs(trailingOnly = TRUE)
tree_file <- args[1]
mode      <- args[2]
selection <- args[3]

cat("================================================\n")
cat("PRUNE TREE\n")
cat("Mode:", mode, "\n")
cat("================================================\n")

tre <- read.tree(tree_file)
stopifnot("Input tree must be rooted" = is.rooted(tre))
was_ultrametric <- !is.null(tre$edge.length) && is.ultrametric(tre, tol = 1e-6)

# Tip labels carry a '.clean' suffix inside the pipeline; accept names with or
# without it so selections can use plain species names.
strip_clean <- function(x) sub("\\.clean$", "", x)
tip_lookup  <- setNames(tre$tip.label, strip_clean(tre$tip.label))

resolve <- function(names_vec) {
  names_vec <- trimws(names_vec)
  names_vec <- names_vec[nzchar(names_vec)]
  hit <- tip_lookup[strip_clean(names_vec)]
  if (any(is.na(hit))) {
    stop("These names are not tips of the tree: ",
         paste(names_vec[is.na(hit)], collapse = ", "),
         "\nTree tips are: ", paste(tre$tip.label, collapse = ", "))
  }
  unname(hit)
}

if (identical(mode, "clade")) {
  seeds <- resolve(strsplit(selection, ",")[[1]])
  if (length(seeds) < 2) {
    stop("--cafe_clade needs at least two tips to define an MRCA.")
  }
  keep <- extract.clade(tre, getMRCA(tre, seeds))$tip.label
  cat("Clade defined by MRCA of:", paste(strip_clean(seeds), collapse = ", "), "\n")
} else if (identical(mode, "species")) {
  wanted <- if (file.exists(selection)) {
    readLines(selection)
  } else {
    strsplit(selection, ",")[[1]]
  }
  wanted <- trimws(wanted)
  keep <- resolve(wanted[nzchar(wanted)])
} else {
  stop("Unknown selection mode: ", mode)
}

keep <- unique(keep)
if (length(keep) < 3) {
  stop("Only ", length(keep), " species selected; CAFE5 needs at least three.")
}
if (length(keep) == Ntip(tre)) {
  cat("NOTE: the selection covers every tip, so the tree is unchanged.\n")
}

pruned <- keep.tip(tre, keep)

stopifnot("Pruned tree must stay rooted" = is.rooted(pruned))
stopifnot("Pruned tree must stay binary" = is.binary(pruned))
if (!is.null(pruned$edge.length)) {
  if (any(!is.finite(pruned$edge.length)) || any(pruned$edge.length <= 0)) {
    stop("Pruned tree contains non-positive or non-finite branch lengths.")
  }
  # Pruning cannot break ultrametricity; check anyway so a calibrated tree is
  # never silently handed to CAFE5 in a state it would reject.
  if (was_ultrametric && !is.ultrametric(pruned, tol = 1e-6)) {
    stop("The input tree was ultrametric but the pruned tree is not.")
  }
}

write.tree(pruned, "SpeciesTree_pruned.nwk")

dropped <- setdiff(tre$tip.label, pruned$tip.label)
write.table(
  data.frame(species = c(sort(pruned$tip.label), sort(dropped)),
             retained = c(rep(TRUE, Ntip(pruned)), rep(FALSE, length(dropped))),
             stringsAsFactors = FALSE),
  "pruned_tree_species.tsv", sep = "\t", quote = FALSE, row.names = FALSE
)

cat("------------------------------------------------\n")
cat("Retained:", Ntip(pruned), "of", Ntip(tre), "species\n")
cat("Dropped: ", length(dropped), "\n")
if (was_ultrametric) {
  d <- node.depth.edgelength(pruned)
  cat("Root age:", round(max(d[seq_len(Ntip(pruned))]), 4), "(tree remains ultrametric)\n")
}
cat("Written: SpeciesTree_pruned.nwk\n")
cat("================================================\n")
