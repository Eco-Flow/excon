#!/usr/bin/env Rscript

# Reference results for tests/go_tables/main.nf.test, computed directly with topGO
# rather than through excon's ChopGO scripts: classic Fisher, every scored term, and
# the number of terms tested in each ontology.
#
# Usage: independent_topgo.R <gff> <gene-to-GO file> <CAFE target list> <CAFE background list>

suppressWarnings(suppressMessages(library(topGO)))

args       <- commandArgs(trailingOnly = TRUE)
gff        <- read.delim(args[1], header = FALSE, comment.char = "#", quote = "")
go         <- read.delim(args[2], header = FALSE, col.names = c("gene", "go"), quote = "")
target     <- readLines(args[3])
background <- readLines(args[4])

genes   <- gff[gff$V3 == "gene", ]
gene_sc <- setNames(genes$V1, sub(";.*", "", sub("^ID=", "", genes$V9)))

test_terms <- function(test, group, selected, universe) {
  ann <- go[go$gene %in% universe, ]
  g2g <- lapply(split(ann$go, ann$gene), unique)
  in_genes <- factor(as.integer(names(g2g) %in% selected), levels = c(0, 1))
  names(in_genes) <- names(g2g)
  do.call(rbind, lapply(c("BP", "MF", "CC"), function(ont) {
    d  <- new("topGOdata", ontology = ont, allGenes = in_genes,
              annot = annFUN.gene2GO, gene2GO = g2g)
    p  <- score(runTest(d, algorithm = "classic", statistic = "fisher"))
    ts <- termStat(d, names(p))
    data.frame(test = test, group = group, ontology = ont, GO.ID = names(p),
               Annotated = ts$Annotated, Significant = ts$Significant,
               P = p, n_tested = length(p))
  }))
}

res <- rbind(
  # Chromosome GO: each scaffold's genes against every GO-annotated gene
  do.call(rbind, lapply(unique(gene_sc), function(sc)
    test_terms("chromo", sc, names(gene_sc)[gene_sc == sc], go$gene))),
  # CAFE GO: the target list against the background list
  test_terms("cafe", "target", target, background)
)

write.table(res, "independent.tsv", sep = "\t", quote = FALSE, row.names = FALSE)
