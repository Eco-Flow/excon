#!/usr/bin/env Rscript

suppressPackageStartupMessages(library(ape))

pdf_dev <- if (capabilities("cairo")) cairo_pdf else pdf
if (!capabilities("cairo")) {
  warning("cairo not available — falling back to pdf(); fonts may not render ",
          "correctly in Inkscape. Add libcairo2-dev to the container to fix this.")
}

args <- commandArgs(trailingOnly = TRUE)
cafe_dir <- if (length(args) >= 1) args[1] else "."

# CAFE writes *_asr.tre as a NEXUS file with per-family trees whose labels carry
# CAFE decorations, e.g. tip "SpA<1>_3" and internal node "<8>*_3" (node id 8).
# read.tree() mis-parses this; parse the first tree and strip the decorations so
# tips are species names and node.label holds the CAFE internal-node ids.
read_cafe_asr <- function(path) {
  lines <- readLines(path, warn = FALSE)
  tl <- grep("^\\s*TREE\\b", lines, ignore.case = TRUE, value = TRUE)
  nwk <- if (length(tl) > 0) sub("^[^=]*=\\s*", "", tl[1]) else paste(lines, collapse = "")
  nwk <- trimws(nwk)
  nwk <- gsub("\\*", "", nwk)
  nwk <- gsub("([A-Za-z0-9_.]+)<[0-9]+>(_-?[0-9]+)?", "\\1", nwk)  # tips
  nwk <- gsub("\\)<([0-9]+)>(_-?[0-9]+)?", ")\\1", nwk)            # internal
  read.tree(text = nwk)
}

tree_files <- list.files(cafe_dir, pattern = "_asr\\.tre$", full.names = TRUE)

if (length(tree_files) == 0) {
  pdf_dev("cafe_node_label_guide.pdf", width = 8, height = 5)
  plot.new()
  text(0.5, 0.5, "No CAFE ASR tree found\n(model may not have converged)",
       cex = 1.2)
  dev.off()
  cat("No *_asr.tre found — blank placeholder written\n")
  quit(status = 0)
}

tree_file <- tree_files[1]
cat("Reading:", tree_file, "\n")
tree <- read_cafe_asr(tree_file)

cat("Tips:", Ntip(tree), "  Internal nodes:", Nnode(tree), "\n")
cat("Node labels:", paste(tree$node.label, collapse = ", "), "\n")

plot_h <- max(6, Ntip(tree) * 0.4 + 2)

pdf_dev("cafe_node_label_guide.pdf", width = 12, height = plot_h)
par(mar = c(2, 1, 4, 10))
plot(tree,
     show.tip.label = TRUE,
     cex            = 0.9,
     label.offset   = 0.002,
     main           = paste0(
       "CAFE Node Number Guide\n",
       "Internal node numbers match CAFE_summary.txt and branch probability tables"
     ))
nodelabels(tree$node.label,
           cex   = 0.75,
           bg    = "#AED6F1",
           col   = "black",
           frame = "rect")
dev.off()

if (requireNamespace("svglite", quietly = TRUE)) {
  svglite::svglite("cafe_node_label_guide.svg", width = 12, height = plot_h)
  par(mar = c(2, 1, 4, 10))
  plot(tree,
       show.tip.label = TRUE,
       cex            = 0.9,
       label.offset   = 0.002,
       main           = paste0(
         "CAFE Node Number Guide\n",
         "Internal node numbers match CAFE_summary.txt and branch probability tables"
       ))
  nodelabels(tree$node.label,
             cex   = 0.75,
             bg    = "#AED6F1",
             col   = "black",
             frame = "rect")
  dev.off()
  cat("Saved: cafe_node_label_guide.svg\n")
} else {
  cat("svglite not available — SVG output skipped\n")
}

cat("Saved: cafe_node_label_guide.pdf\n")
