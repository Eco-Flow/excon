#!/usr/bin/env Rscript

# Time-calibrate the species tree with IQ-TREE's least-square dating (LSD2),
# --dating_method lsd2. The calibration table is the one --tree_calibrations
# takes for ape::chronos (clade, tips, age_min, age_max; see date_tree.R).
#
#   lsd2_dating.R prepare <rooted tree> <calibrations.tsv>
#       Checks every calibration against the tree before dating: all tips present,
#       no two calibrations on the same node, 0 < age_min <= age_max, nested ages
#       that can hold, and at least one exact age (LSD2 has no unique solution from
#       bounds alone). Writes
#         lsd2_dates.txt    IQ-TREE --date file: "tip1,tip2<TAB>-age" for an exact age,
#                           "-age_max:-age_min" for bounds (dates, so the past is negative)
#         lsd2_outgroup.txt the tips on the smaller side of the root, for -o, so LSD2
#                           keeps the tree's root rather than estimating a new one
#
#   lsd2_dating.R check <dated.timetree.nex> <calibrations.tsv> <rooted tree> <iqtree log>
#       Checks the dated tree: LSD2 reported no warning, the tree is ultrametric and
#       binary with positive branch lengths, has the input's topology and root, and
#       every calibrated node lies within its bounds. Writes SpeciesTree_dated.nwk,
#       dating_calibrations.tsv and dating_qc.tsv (as date_tree.R does).

suppressPackageStartupMessages(library(ape))

args <- commandArgs(trailingOnly = TRUE)
mode <- args[1]

strip_clean <- function(x) sub("\\.clean$", "", x)

read_calibrations <- function(path) {
  calib <- read.delim(path, stringsAsFactors = FALSE)
  missing_cols <- setdiff(c("clade", "tips", "age_min", "age_max"), names(calib))
  if (length(missing_cols)) stop("Calibration table is missing column(s): ", paste(missing_cols, collapse = ", "))
  if (!nrow(calib)) stop("Calibration table is empty.")
  calib
}

# Tip labels may carry the pipeline's '.clean' suffix; calibrations may omit it
resolve_tips <- function(tre, spec, clade) {
  lookup <- setNames(tre$tip.label, strip_clean(tre$tip.label))
  wanted <- trimws(strsplit(spec, ",")[[1]])
  wanted <- wanted[nzchar(wanted)]
  if (length(wanted) < 2) stop("Calibration '", clade, "' needs at least two tips to define an MRCA.")
  resolved <- lookup[strip_clean(wanted)]
  if (any(is.na(resolved))) {
    stop("Calibration '", clade, "' names tips absent from the tree: ",
         paste(wanted[is.na(resolved)], collapse = ", "),
         "\nTree tips are: ", paste(tre$tip.label, collapse = ", "))
  }
  unname(resolved)
}

calibrated_nodes <- function(tre, calib) {
  vapply(seq_len(nrow(calib)), function(i)
    getMRCA(tre, resolve_tips(tre, calib$tips[i], calib$clade[i])), integer(1))
}

node_age <- function(tree, node) {
  d <- node.depth.edgelength(tree)
  max(d[seq_len(Ntip(tree))]) - d[node]
}

if (mode == "prepare") {
  tre   <- read.tree(args[2])
  calib <- read_calibrations(args[3])
  if (!is.rooted(tre)) stop("The species tree must be rooted.")
  tre$node.label <- NULL

  nodes <- calibrated_nodes(tre, calib)
  if (anyDuplicated(nodes)) {
    dup <- calib$clade[duplicated(nodes) | duplicated(nodes, fromLast = TRUE)]
    stop("These calibrations resolve to the same node: ", paste(dup, collapse = ", "))
  }
  bad <- calib$age_min > calib$age_max | calib$age_min <= 0
  if (any(bad)) stop("Calibrations need 0 < age_min <= age_max: ", paste(calib$clade[bad], collapse = ", "))
  tips_under <- lapply(nodes, function(n) extract.clade(tre, n)$tip.label)
  for (i in seq_along(nodes)) for (j in seq_along(nodes)) {
    if (i != j && all(tips_under[[j]] %in% tips_under[[i]]) && calib$age_max[i] <= calib$age_min[j]) {
      stop(sprintf("Calibrations '%s' (%.2f-%.2f Ma) and '%s' (%.2f-%.2f Ma) conflict: '%s' lies inside '%s', so it must be younger, but its minimum age is not below the maximum age of '%s'.",
                   calib$clade[i], calib$age_min[i], calib$age_max[i],
                   calib$clade[j], calib$age_min[j], calib$age_max[j],
                   calib$clade[j], calib$clade[i], calib$clade[i]))
    }
  }
  if (!any(calib$age_min == calib$age_max)) {
    stop("LSD2 needs at least one calibration with an exact age (age_min equal to age_max); ",
         "from bounds alone it has no unique solution. Fix one age, or use --dating_method chronos.")
  }

  # IQ-TREE date file: dates, not ages, with the tips at 0 (--date-tip 0)
  fmt <- function(x) format(x, scientific = FALSE, trim = TRUE)
  dates <- ifelse(calib$age_min == calib$age_max,
                  paste0("-", fmt(calib$age_min)),
                  paste0("-", fmt(calib$age_max), ":-", fmt(calib$age_min)))
  # The full tip list of each node, so the MRCA is exactly the calibrated node
  lines <- vapply(seq_along(nodes), function(i)
    paste0(paste(tips_under[[i]], collapse = ","), "\t", dates[i]), character(1))
  writeLines(lines, "lsd2_dates.txt")

  # Keep the input root: name the tips on the smaller side of it as the outgroup
  root  <- Ntip(tre) + 1
  sides <- tre$edge[tre$edge[, 1] == root, 2]
  side_tips <- lapply(sides, function(n) if (n <= Ntip(tre)) tre$tip.label[n] else extract.clade(tre, n)$tip.label)
  writeLines(paste(side_tips[[which.min(lengths(side_tips))]], collapse = ","), "lsd2_outgroup.txt")

  cat("Calibrations checked against the tree:\n")
  for (i in seq_along(nodes)) {
    cat(sprintf("   %-28s %d tips   %s Ma\n", calib$clade[i], length(tips_under[[i]]),
                if (calib$age_min[i] == calib$age_max[i]) fmt(calib$age_min[i])
                else paste0(fmt(calib$age_min[i]), "-", fmt(calib$age_max[i]))))
  }

} else if (mode == "check") {
  calib <- read_calibrations(args[3])
  input <- read.tree(args[4])
  log   <- readLines(args[5], warn = FALSE)

  # LSD2 reports a non-unique or failed solution as a warning; treat it as fatal,
  # as date_tree.R does for chronos
  # (LSD2's own warnings only; IQ-TREE logs unrelated "WARNING:" lines too)
  warn <- grep("^\\*WARNINGS|not enough input date information", log, value = TRUE)
  if (length(warn)) stop("LSD2 reported: ", paste(unique(trimws(warn)), collapse = " | "))

  # IQ-TREE's time tree in NEXUS, branch lengths in time units, nodes annotated with
  # [&date="..."]; drop the annotations to read it as plain Newick
  nex <- readLines(args[2], warn = FALSE)
  tl  <- grep("^\\s*tree\\s", nex, ignore.case = TRUE, value = TRUE)
  if (!length(tl)) stop("No tree found in ", args[2])
  nwk <- gsub("\\[&[^]]*\\]", "", sub("^[^=]*=\\s*", "", tl[1]))
  fit <- read.tree(text = nwk)

  if (!is.ultrametric(fit, tol = 1e-6)) stop("The dated tree is not ultrametric.")
  if (!is.binary(fit)) stop("The dated tree is not binary.")
  if (any(!is.finite(fit$edge.length)) || any(fit$edge.length <= 0)) {
    stop("The dated tree has non-positive or non-finite branch lengths: ",
         "a calibration probably forces two nodes to the same age.")
  }
  input$node.label <- NULL
  if (!setequal(fit$tip.label, input$tip.label)) stop("The dated tree's tips differ from the input tree's.")
  if (!all.equal.phylo(fit, input, use.edge.length = FALSE, use.tip.label = TRUE)) {
    stop("LSD2 changed the tree's topology or root.")
  }

  nodes <- calibrated_nodes(fit, calib)
  ages  <- vapply(nodes, function(n) node_age(fit, n), numeric(1))
  tol    <- 1e-3
  within <- ages >= calib$age_min - tol & ages <= calib$age_max + tol
  if (!all(within)) {
    stop("The dated tree does not reproduce these calibrations: ",
         paste(sprintf("%s (wanted %.2f-%.2f, got %.2f)", calib$clade[!within],
                       calib$age_min[!within], calib$age_max[!within], ages[!within]),
               collapse = "; "))
  }

  write.tree(fit, "SpeciesTree_dated.nwk")
  write.table(data.frame(clade = calib$clade, tips = calib$tips, node = nodes,
                         age_min_Ma = calib$age_min, age_max_Ma = calib$age_max,
                         fitted_age_Ma = round(ages, 6), within_bounds = within),
              "dating_calibrations.tsv", sep = "\t", quote = FALSE, row.names = FALSE)
  rate <- sub(".*rate ([^,]+),.*", "\\1", grep("rate .*tMRCA", log, value = TRUE)[1])
  depths <- node.depth.edgelength(fit)
  write.table(data.frame(method = "lsd2 (IQ-TREE --date)", substitution_rate = rate,
                         n_tips = Ntip(fit), n_calibrations = nrow(calib),
                         root_age_Ma = round(max(depths[seq_len(Ntip(fit))]), 6),
                         min_edge_Ma = round(min(fit$edge.length), 8),
                         ultrametric = TRUE, binary = TRUE),
              "dating_qc.tsv", sep = "\t", quote = FALSE, row.names = FALSE)
  cat("Dated tree written: SpeciesTree_dated.nwk (root age", round(max(depths), 3), "Ma)\n")

} else {
  stop("Usage: lsd2_dating.R prepare <tree> <calibrations> | check <timetree.nex> <calibrations> <tree> <iqtree log>")
}
