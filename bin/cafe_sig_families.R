#!/usr/bin/env Rscript

# Per-node significant-family tables + retrieval of alignments/gene trees for the
# families significantly expanded at the internal node(s) of biological interest.
#
# Reads the best CAFE5 output directory (Base_*_branch_probabilities.tab,
# Base_*_change.tab and Base_*_asr.tre) and the OrthoFinder N0.tsv, and writes:
#   * changes_per_node.tsv             — every branch with a non-zero change
#   * significant_changes_per_node.tsv — the subset with branch p-value <= cutoff
#   * (if focus clades + OrthoFinder dirs supplied) focus_families/<node>/... plus
#     focus_families_manifest.tsv, holding the MSA (.fa) and gene tree of every
#     family SIGNIFICANTLY EXPANDED at each focus node.
#
# Positional args (from cafe_sig_families.nf):
#   args[1] = CAFE results directory (e.g. Out_cafe)
#   args[2] = N0.tsv path
#   args[3] = p-value cutoff (default 0.05)
#   args[4] = focus clade spec, '|'-separated. Each clade is either a bare CAFE
#             node label (integer, as in cafe_node_label_guide.pdf) or a
#             comma-separated set of >=2 tip species whose MRCA defines the node.
#             "" / "NA" to skip focus extraction.
#   args[5] = OrthoFinder MSA directory (per-OG *.fa), or "NA"/missing to skip
#   args[6] = OrthoFinder gene-tree directory (per-OG *_tree.txt), or "NA"/missing

suppressPackageStartupMessages({
  library(ape)
  library(data.table)
})

args      <- commandArgs(trailingOnly = TRUE)
cafe_dir  <- if (length(args) >= 1) args[1] else "."
n0_path   <- if (length(args) >= 2) args[2] else "N0.tsv"
pcut      <- if (length(args) >= 3) as.numeric(args[3]) else 0.05
focus_arg <- if (length(args) >= 4) args[4] else ""
msa_dir   <- if (length(args) >= 5) args[5] else "NA"
tree_dir  <- if (length(args) >= 6) args[6] else "NA"

is_na <- function(x) is.null(x) || is.na(x) || x %in% c("", "NA", "NO_FILE")

# CAFE writes *_asr.tre as a NEXUS file with per-family trees whose labels carry
# CAFE decorations, e.g. tip "SpA<1>_3" and internal node "<8>*_3" (node id 8).
# Return a plain phylo whose tips are the species names and whose node.label
# entries are the CAFE internal-node ids used in the .tab column headers.
read_cafe_asr <- function(path) {
  lines <- readLines(path, warn = FALSE)
  tl <- grep("^\\s*TREE\\b", lines, ignore.case = TRUE, value = TRUE)
  nwk <- if (length(tl) > 0) sub("^[^=]*=\\s*", "", tl[1]) else paste(lines, collapse = "")
  nwk <- trimws(nwk)
  nwk <- gsub("\\*", "", nwk)                                      # significance stars
  nwk <- gsub("([A-Za-z0-9_.]+)<[0-9]+>(_-?[0-9]+)?", "\\1", nwk)  # tips  -> species name
  nwk <- gsub("\\)<([0-9]+)>(_-?[0-9]+)?", ")\\1", nwk)            # internal -> node id
  read.tree(text = nwk)
}

## ------------------------------------------------------------------
## Locate the CAFE output tables
## ------------------------------------------------------------------
f_prob <- list.files(cafe_dir, pattern = "_branch_probabilities\\.tab$", full.names = TRUE)
f_chg  <- list.files(cafe_dir, pattern = "_change\\.tab$",               full.names = TRUE)
f_asr  <- list.files(cafe_dir, pattern = "_asr\\.tre$",                  full.names = TRUE)

if (length(f_prob) == 0 || length(f_chg) == 0) {
  # Model did not converge / no tables — write empty outputs so the pipeline continues.
  fwrite(data.table(Node = character(), HOG = character(), change = integer(),
                    pvalue = numeric(), direction = character(), significant = logical()),
         "changes_per_node.tsv", sep = "\t")
  fwrite(data.table(Node = character(), HOG = character(), change = integer(),
                    pvalue = numeric(), direction = character()),
         "significant_changes_per_node.tsv", sep = "\t")
  cat("No CAFE branch tables found in", cafe_dir, "— empty outputs written.\n")
  quit(status = 0)
}

prob <- fread(f_prob[1], header = TRUE, check.names = FALSE)
chg  <- fread(f_chg[1],  header = TRUE, check.names = FALSE)

# First column is the family id (#FamilyID / FamilyID)
id_name <- names(prob)[1]
setnames(prob, id_name, "HOG")
setnames(chg,  names(chg)[1], "HOG")

# Map each CAFE column header to a node id: internal nodes -> "<24>" become "24",
# tips -> "Species<1>" become "Species".
node_id_of <- function(col) {
  if (grepl("^<[0-9]+>$", col)) return(sub("^<([0-9]+)>$", "\\1", col))
  if (grepl("<[0-9]+>$", col))  return(sub("<[0-9]+>$", "", col))
  col
}
node_cols <- setdiff(names(prob), "HOG")
node_ids  <- vapply(node_cols, node_id_of, character(1))

parse_change <- function(x) {
  x <- gsub("[+]", "", x)
  suppressWarnings(as.integer(x))
}

## ------------------------------------------------------------------
## Long tables: one row per (family, node)
## ------------------------------------------------------------------
prob_l <- melt(prob, id.vars = "HOG", variable.name = "col", value.name = "pvalue")
chg_l  <- melt(chg,  id.vars = "HOG", variable.name = "col", value.name = "change_raw")
prob_l[, Node := node_ids[as.character(col)]]
chg_l[,  Node := node_ids[as.character(col)]]
prob_l[, pvalue := suppressWarnings(as.numeric(pvalue))]
chg_l[,  change := parse_change(change_raw)]

long <- merge(chg_l[, .(HOG, Node, change)],
              prob_l[, .(HOG, Node, pvalue)],
              by = c("HOG", "Node"), all = TRUE)
long <- long[!is.na(change) & change != 0]
long[, direction   := ifelse(change > 0, "expansion", "contraction")]
long[, significant := !is.na(pvalue) & pvalue <= pcut]

setorder(long, Node, pvalue, -change)
fwrite(long[, .(Node, HOG, change, pvalue, direction, significant)],
       "changes_per_node.tsv", sep = "\t")
fwrite(long[significant == TRUE, .(Node, HOG, change, pvalue, direction)],
       "significant_changes_per_node.tsv", sep = "\t")
cat("Wrote changes_per_node.tsv (", nrow(long), " changed branches) and ",
    "significant_changes_per_node.tsv (", sum(long$significant), " significant).\n", sep = "")

## ------------------------------------------------------------------
## Focus-node extraction of alignments / gene trees
## ------------------------------------------------------------------
if (is_na(focus_arg)) {
  cat("No --cafe_focus_clades supplied — skipping alignment/gene-tree retrieval.\n")
  quit(status = 0)
}

# Resolve each focus clade spec to a CAFE node id
resolve_focus_nodes <- function(spec, asr_file) {
  clades <- strsplit(spec, "|", fixed = TRUE)[[1]]
  clades <- trimws(clades[clades != ""])
  res <- list()
  tree <- NULL
  for (cl in clades) {
    if (grepl("^[0-9]+$", cl)) {
      res[[length(res) + 1]] <- list(node = cl, label = paste0("node", cl))
      next
    }
    sp <- trimws(strsplit(cl, ",")[[1]])
    if (length(sp) < 2) {
      warning("Focus clade '", cl, "' needs a node number or >=2 species — skipped.")
      next
    }
    if (is.null(tree)) {
      if (length(asr_file) == 0) {
        warning("No asr tree to resolve species MRCA.")
        return(res)
      }
      tree <- read_cafe_asr(asr_file[1])
    }
    idx <- match(sp, tree$tip.label)
    if (any(is.na(idx))) {
      warning("Focus species not found in tree: ",
              paste(sp[is.na(idx)], collapse = ", "), " — clade skipped.")
      next
    }
    mrca <- getMRCA(tree, tree$tip.label[idx])
    lbl  <- tree$node.label[mrca - Ntip(tree)]
    node <- gsub("[^0-9]", "", lbl)
    res[[length(res) + 1]] <- list(node = node,
                                   label = paste0("node", node, "_",
                                                  paste(sp, collapse = "-")))
  }
  res
}

focus <- resolve_focus_nodes(focus_arg, f_asr)
if (length(focus) == 0) {
  cat("No focus nodes could be resolved — skipping alignment retrieval.\n")
  quit(status = 0)
}

# HOG -> OG map (OrthoFinder v2 has an OG column; v3 the ids are already OGs)
n0 <- fread(n0_path, header = TRUE, check.names = FALSE)
hog_col <- if ("HOG" %in% names(n0)) "HOG" else "Orthogroup"
if ("OG" %in% names(n0)) {
  hog2og <- setNames(n0[["OG"]], n0[[hog_col]])
} else {
  hog2og <- setNames(n0[[hog_col]], n0[[hog_col]])
}

dir.create("focus_families", showWarnings = FALSE)
manifest <- data.table()

have_msa  <- !is_na(msa_dir)  && dir.exists(msa_dir)
have_tree <- !is_na(tree_dir) && dir.exists(tree_dir)
if (!have_msa)  cat("MSA directory not available — alignments will not be copied.\n")
if (!have_tree) cat("Gene-tree directory not available — trees will not be copied.\n")

find_one <- function(dir, og) {
  # Match OrthoFinder file for an OG id, tolerant of naming (OG.., *_tree.txt, .fa)
  hits <- list.files(dir, pattern = paste0("(^|[^0-9])", og, "([^0-9]|$|_|\\.)"),
                     full.names = TRUE)
  if (length(hits) == 0) return(NA_character_) else hits[1]
}

for (fc in focus) {
  node <- fc$node
  sub  <- long[Node == node & direction == "expansion" & significant == TRUE]
  if (nrow(sub) == 0) {
    cat("Focus node", node, "(", fc$label, "): no significantly expanded families.\n")
    next
  }
  outdir <- file.path("focus_families", fc$label)
  dir.create(outdir, showWarnings = FALSE, recursive = TRUE)
  for (i in seq_len(nrow(sub))) {
    hog <- sub$HOG[i]
    og  <- hog2og[[hog]]
    if (is.null(og) || is.na(og)) og <- hog
    msa_src  <- if (have_msa)  find_one(msa_dir, og)  else NA_character_
    tree_src <- if (have_tree) find_one(tree_dir, og) else NA_character_
    if (!is.na(msa_src))  file.copy(msa_src,  file.path(outdir, basename(msa_src)),  overwrite = TRUE)
    if (!is.na(tree_src)) file.copy(tree_src, file.path(outdir, basename(tree_src)), overwrite = TRUE)
    manifest <- rbind(manifest, data.table(
      focus_node = node, focus_label = fc$label, HOG = hog, OG = og,
      change = sub$change[i], pvalue = sub$pvalue[i],
      msa_file  = if (is.na(msa_src))  "" else basename(msa_src),
      tree_file = if (is.na(tree_src)) "" else basename(tree_src),
      msa_found  = !is.na(msa_src),
      tree_found = !is.na(tree_src)
    ))
  }
  cat("Focus node", node, "(", fc$label, "): ", nrow(sub),
      " significantly expanded families extracted.\n", sep = "")
}

if (nrow(manifest) > 0) {
  fwrite(manifest, "focus_families_manifest.tsv", sep = "\t")
  cat("Wrote focus_families_manifest.tsv (", nrow(manifest), " families).\n", sep = "")
} else {
  fwrite(data.table(focus_node = character(), focus_label = character(),
                    HOG = character(), OG = character(), change = integer(),
                    pvalue = numeric(), msa_file = character(),
                    tree_file = character(), msa_found = logical(),
                    tree_found = logical()),
         "focus_families_manifest.tsv", sep = "\t")
  cat("No focus families extracted — empty manifest written.\n")
}
