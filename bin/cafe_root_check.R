#!/usr/bin/env Rscript
# cafe_root_check.R
# Warns when one species sits alone on one side of the root of the CAFE5 tree.
#
# CAFE5 infers the root family size, and with a single species on one side, a
# change on that species' branch, or on the branch leading to all the others,
# can't be told apart from a different root size. The counts on those two
# branches then follow the root-size prior (uniform vs Poisson) more than the
# data: e.g. with a Poisson root, both can show expansions only.
#
# Usage: cafe_root_check.R <CAFE5 results dir>
# Writes root_split_warning.txt (and prints it) only when the root is split that way.

suppressPackageStartupMessages(library(ape))

cafe_dir   <- commandArgs(trailingOnly = TRUE)[1]
tree_files <- list.files(cafe_dir, pattern = "_asr\\.tre$", full.names = TRUE)
if (length(tree_files) == 0) quit(status = 0)

# *_asr.tre is NEXUS, one tree per family, all with the species topology
tree <- read.nexus(tree_files[1])
if (inherits(tree, "multiPhylo")) tree <- tree[[1]]

# CAFE5 labels: tips "Species<1>*_3", internal nodes "<14>*_8"
cafe_id  <- function(label) sub("^.*<([0-9]+)>.*$", "\\1", label)
species  <- sub("<.*$", "", tree$tip.label)
n_tip    <- Ntip(tree)
root     <- n_tip + 1
children <- tree$edge[tree$edge[, 1] == root, 2]
lone     <- children[children <= n_tip]

if (length(children) == 2 && length(lone) == 1) {
  other <- setdiff(children, lone)
  other_label <- if (other <= n_tip) tree$tip.label[other] else tree$node.label[other - n_tip]
  msg <- c(
    sprintf("WARNING: %s is alone on one side of the root of the CAFE5 tree.", species[lone]),
    sprintf("Changes on its branch (CAFE5 node %s) and on the branch to the other %d species (node %s)",
            cafe_id(tree$tip.label[lone]), n_tip - 1, cafe_id(other_label)),
    "can't be told apart from a different root family size, so the expansions and contractions",
    "reported on these two branches mostly reflect the root-size prior (uniform vs Poisson) rather",
    "than evidence. Interpret them with caution; another species on that side of the root fixes this."
  )
  writeLines(msg, "root_split_warning.txt")
  message(paste(msg, collapse = "\n"))
}
