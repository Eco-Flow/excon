#!/usr/bin/env Rscript

suppressPackageStartupMessages(library(ape))

args <- commandArgs(trailingOnly = TRUE)
cafe_dir <- if (length(args) >= 1) args[1] else "."

tree_files <- list.files(cafe_dir, pattern = "_asr\\.tre$", full.names = TRUE)

if (length(tree_files) == 0) {
    pdf("cafe_node_label_guide.pdf", width = 8, height = 5)
    plot.new()
    text(0.5, 0.5, "No CAFE ASR tree found\n(model may not have converged)", cex = 1.2)
    dev.off()
    cat("No *_asr.tre found — blank placeholder written\n")
    quit(status = 0)
}

tree_file <- tree_files[1]
cat("Reading:", tree_file, "\n")
tree <- read.tree(tree_file)

cat("Tips:", Ntip(tree), "  Internal nodes:", Nnode(tree), "\n")
cat("Node labels:", paste(tree$node.label, collapse = ", "), "\n")

plot_h <- max(6, Ntip(tree) * 0.4 + 2)

pdf("cafe_node_label_guide.pdf", width = 12, height = plot_h)
par(mar = c(2, 1, 4, 10))
plot(tree,
     show.tip.label = TRUE,
     cex            = 0.9,
     label.offset   = 0.002,
     main           = "CAFE Node Number Guide\nInternal node numbers match CAFE_summary.txt and branch probability tables")
nodelabels(tree$node.label,
           cex   = 0.75,
           bg    = "#AED6F1",
           col   = "black",
           frame = "rect")
dev.off()

cat("Saved: cafe_node_label_guide.pdf\n")
