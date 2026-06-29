#!/usr/bin/env Rscript
# cafe_plot_altviz.R
# Two alternative CAFE result tree figures derived from CAFE_summary.txt:
#   1. cafe_sig_hog_tree  — significantly evolving gene families per branch (p <= 0.05)
#   2. cafe_sig_gene_tree — net genes gained/lost in those significant families
#
# The standard cafeplotter output shows ALL families (no p-value filter).
# These figures show only the statistically significant signal.
#
# Usage: cafe_plot_altviz.R <cafe_dir> <CAFE_summary.txt>

suppressPackageStartupMessages(library(ape))

args         <- commandArgs(trailingOnly = TRUE)
cafe_dir     <- if (length(args) >= 1) args[1] else "."
summary_file <- if (length(args) >= 2) args[2] else "CAFE_summary.txt"

# ── Tree ───────────────────────────────────────────────────────────────────────
tree_files <- list.files(cafe_dir, pattern = "_asr\\.tre$", full.names = TRUE)
if (length(tree_files) == 0) {
  message("No *_asr.tre found — skipping cafe_plot_altviz")
  quit(status = 0)
}
tree <- read.tree(tree_files[1])
# CAFE5 *_asr.tre holds one reconstructed tree per gene family (a multiPhylo).
# They all share the species topology, so use the first for the tree scaffold.
if (inherits(tree, "multiPhylo")) {
  message("asr file contains ", length(tree), " trees — using the first for topology")
  tree <- tree[[1]]
}
# Strip angle brackets and leading 'm' that CAFE5 sometimes puts on node labels
tree$node.label <- gsub("^[<m]|[<>]$", "", tree$node.label)

# ── CAFE_summary.txt ───────────────────────────────────────────────────────────
if (!file.exists(summary_file)) {
  message("Not found: ", summary_file); quit(status = 0)
}
dat <- read.table(summary_file, header = TRUE, sep = "\t",
                  comment.char = "#", quote = "", stringsAsFactors = FALSE,
                  check.names = FALSE)
colnames(dat)[1] <- "node"

# ── Map CAFE_summary rows to ape node indices ─────────────────────────────────
# Tips:            ape indices 1 .. n_tip      (tree$tip.label)
# Internal nodes:  ape indices n_tip+1 .. N    (tree$node.label[i] = node n_tip+i)
n_tip  <- Ntip(tree)
n_node <- Nnode(tree)
N      <- n_tip + n_node

mk <- function() rep(NA_real_, N)
exp_hog    <- mk()
contr_hog  <- mk()
exp_gene   <- mk()
contr_gene <- mk()

for (i in seq_len(nrow(dat))) {
  nm <- as.character(dat$node[i])
  if (nm %in% tree$tip.label) {
    idx <- which(tree$tip.label == nm)
  } else {
    m <- which(tree$node.label == nm)
    if (!length(m)) next
    idx <- n_tip + m[1]
  }
  exp_hog[idx]    <- dat$Expansion_HOGs_significant[i]
  contr_hog[idx]  <- dat$Contraction_HOGs_significant[i]          # positive count
  exp_gene[idx]   <- dat$Expansion_genes_significant[i]           # positive sum
  contr_gene[idx] <- abs(dat$Contraction_genes_significant[i])    # make positive
}

# ── Plot helper ───────────────────────────────────────────────────────────────
pdf_dev <- if (capabilities("cairo")) cairo_pdf else pdf
if (!capabilities("cairo"))
  warning("cairo not available — falling back to pdf(); fonts may not render well in Inkscape")

plot_h <- max(6, n_tip * 0.45 + 2.5)

do_plot <- function(exp_v, contr_v, title) {
  tree_len <- max(node.depth.edgelength(tree))

  # Tip labels: "Species name   +exp / −contr"
  tip_text <- mapply(function(sp, e, c) {
    if (is.na(e)) sp
    else sprintf("%s   +%g / −%g", sp, e, c)
  }, tree$tip.label, exp_v[seq_len(n_tip)], contr_v[seq_len(n_tip)])

  par(mar = c(2, 1, 3.5, 1))
  plot(tree,
       show.tip.label = FALSE,
       main     = title,
       cex.main = 0.9,
       x.lim    = c(0, tree_len * 1.85))

  tiplabels(text   = tip_text,
            frame  = "none",
            adj    = c(0, 0.5),
            cex    = 0.75,
            offset = tree_len * 0.02)

  # Internal node labels — only for nodes that matched a CAFE_summary row
  int_exp   <- exp_v[seq(n_tip + 1, N)]
  int_contr <- contr_v[seq(n_tip + 1, N)]
  valid_i   <- which(!is.na(int_exp))
  if (length(valid_i) > 0) {
    node_text <- sprintf("+%g / −%g", int_exp[valid_i], int_contr[valid_i])
    nodelabels(text  = node_text,
               node  = valid_i + n_tip,
               frame = "rect",
               bg    = "#EAF4FB",
               col   = "#1A5276",
               cex   = 0.62)
  }
}

save_both <- function(stem, exp_v, contr_v, title) {
  pdf_dev(paste0(stem, ".pdf"), width = 14, height = plot_h)
  do_plot(exp_v, contr_v, title)
  dev.off()

  if (requireNamespace("svglite", quietly = TRUE)) {
    svglite::svglite(paste0(stem, ".svg"), width = 14, height = plot_h)
    do_plot(exp_v, contr_v, title)
    dev.off()
    message("Saved: ", stem, ".pdf / .svg")
  } else {
    message("Saved: ", stem, ".pdf  (svglite not available — SVG skipped)")
  }
}

# Figures are supplementary — never let a plotting error abort the pipeline.
ok <- tryCatch({
  # ── Figure 1: Significant HOGs ──────────────────────────────────────────────
  save_both(
    "cafe_sig_hog_tree",
    exp_hog, contr_hog,
    "Significantly evolving gene families per branch (+expanded / −contracted, p ≤ 0.05)"
  )

  # ── Figure 2: Net genes in significant families ─────────────────────────────
  save_both(
    "cafe_sig_gene_tree",
    exp_gene, contr_gene,
    "Net genes in significant families per branch (+gained / −lost, p ≤ 0.05)"
  )
  TRUE
}, error = function(e) {
  message("cafe_plot_altviz: figure generation failed (", conditionMessage(e),
          ") — skipping (non-fatal)")
  FALSE
})

message(if (ok) "Done." else "Skipped alt-viz figures (non-fatal).")
quit(status = 0)
