#!/usr/bin/env Rscript
# plot_go_enrichment.R
# Professional GO enrichment visualisation: bar plot + dot plot
# Outputs PDF and SVG for each plot type.
#
# Usage:
#   Rscript plot_go_enrichment.R <topgo_result_table> \
#       [min_p] [num_results] [open_plot] [correction] \
#       [enrich_cutoff] [min_gene] [max_gene] [inferred_go_file]
#
# Arguments (all optional beyond the first):
#   min_p          max p-value to display          (default 0.05)
#   num_results    max GO terms per ontology        (default 10)
#   open_plot      open PDF after plotting 0|1      (default 0)
#   correction     p-value column name to use       (default bonferroni)
#   enrich_cutoff  minimum fold-enrichment          (default 0)
#   min_gene       minimum annotated genes           (default 5)
#   max_gene       maximum annotated genes           (default 100000)
#   inferred_go    path to GO inference filter file  (optional)

suppressPackageStartupMessages({
  for (pkg in c("ggplot2", "dplyr", "stringr")) {
    if (!requireNamespace(pkg, quietly = TRUE))
      stop("Required R package missing: ", pkg,
           "\nInstall with: install.packages('", pkg, "')")
    library(pkg, character.only = TRUE)
  }
})

# GO.db for full untruncated term names (topGO truncates at ~40 chars)
has_godb <- requireNamespace("GO.db", quietly = TRUE) &&
            requireNamespace("AnnotationDbi", quietly = TRUE)
if (has_godb) {
  suppressPackageStartupMessages({
    library(GO.db)
    library(AnnotationDbi)
  })
}

# ── Parse arguments ────────────────────────────────────────────────────────────
args <- commandArgs(trailingOnly = TRUE)
if (length(args) < 1)
  stop("Usage: Rscript plot_go_enrichment.R <topgo_result_table> [options...]")

input_file    <- args[1]
min_p         <- if (length(args) >= 2) as.numeric(args[2]) else 0.05
num_results   <- if (length(args) >= 3) as.integer(args[3]) else 10
open_plot     <- if (length(args) >= 4) as.integer(args[4]) else 0
correction    <- if (length(args) >= 5) args[5] else "bonferroni"
enrich_cutoff <- if (length(args) >= 6) as.numeric(args[6]) else 0
min_gene      <- if (length(args) >= 7) as.integer(args[7]) else 5
max_gene      <- if (length(args) >= 8) as.integer(args[8]) else 100000
inferred_file <- if (length(args) >= 9 && nchar(args[9]) > 0) args[9] else NULL

# ── Load data ──────────────────────────────────────────────────────────────────
if (!file.exists(input_file)) stop("Input file not found: ", input_file)

dat <- read.table(input_file, header = TRUE, sep = "\t", quote = "",
                  fill = TRUE, stringsAsFactors = FALSE, comment.char = "")

# ── Inferred GO filter ─────────────────────────────────────────────────────────
# Keeps only terms present in the inference file (column 2 = GO ID)
if (!is.null(inferred_file) && file.exists(inferred_file)) {
  inf      <- read.table(inferred_file, sep = "\t", header = FALSE,
                         stringsAsFactors = FALSE, fill = TRUE)
  keep_ids <- inf[[2]]
  dat      <- dat[dat$GO.ID %in% keep_ids, ]
}

# ── Replace truncated topGO names with full names from GO.db ──────────────────
if (has_godb) {
  go_ids   <- dat$GO.ID
  go_terms <- suppressMessages(
    AnnotationDbi::select(GO.db, keys = go_ids,
                          columns = "TERM", keytype = "GOID")
  )
  # match back by position (select returns 1:1 when keys are valid GO IDs)
  m <- match(dat$GO.ID, go_terms$GOID)
  full_terms <- go_terms$TERM[m]
  dat$Term <- ifelse(!is.na(full_terms), full_terms, dat$Term)
}

# ── Sanitise ───────────────────────────────────────────────────────────────────
dat$Term <- ifelse(is.na(dat$Term) | dat$Term == "", dat$GO.ID, dat$Term)
dat$Term <- gsub("'", "", dat$Term)

# ── Append GO ID to label: "full term name (GO:XXXXXXX)" ──────────────────────
dat$Term <- paste0(dat$Term, "\n(", dat$GO.ID, ")")

dat$FoldChange <- suppressWarnings(as.numeric(
  gsub("inf", "100", tolower(as.character(dat$FoldChange)))
))
dat$FoldChange[is.infinite(dat$FoldChange) | is.na(dat$FoldChange)] <- 100

# Resolve chosen p-value column; handle "< 1e-30" strings
if (!correction %in% colnames(dat))
  stop("Correction '", correction, "' not found. Available columns: ",
       paste(colnames(dat), collapse = ", "))

raw_p       <- gsub("< *1e-30", "1e-30", as.character(dat[[correction]]))
dat$pval    <- suppressWarnings(as.numeric(raw_p))
dat$pval[is.na(dat$pval) | dat$pval == 0] <- 1e-30

# ── Filter ─────────────────────────────────────────────────────────────────────
dat <- dat[
  dat$Annotated > min_gene  &
  dat$Annotated < max_gene  &
  dat$FoldChange >= enrich_cutoff &
  dat$pval < min_p,
]

if (nrow(dat) == 0) {
  message("No significant GO terms with these thresholds. No plots produced.")
  quit(status = 0)
}

# ── Derived columns ────────────────────────────────────────────────────────────
dat$logP <- -log10(dat$pval)

# ── Top N per ontology ─────────────────────────────────────────────────────────
valid_ont <- c("BP", "MF", "CC")
dat <- dat %>%
  filter(ontology %in% valid_ont) %>%
  group_by(ontology) %>%
  arrange(pval) %>%
  slice_head(n = num_results) %>%
  ungroup() %>%
  mutate(ontology = factor(ontology, levels = valid_ont))

if (nrow(dat) == 0) {
  message("No terms remain after ontology filtering.")
  quit(status = 0)
}

# ── Order terms for display: within each ontology, ascending logP ─────────────
# Wrapping is applied at the scale layer (scale_y_discrete) — NOT on the data.
# This matches enrichplot's approach and avoids broken factor levels with \n.
dat <- dat %>%
  arrange(ontology, logP) %>%
  mutate(Term_label = factor(Term, levels = unique(Term)))

# Label wrapping closure (enrichplot pattern): 30-char word-boundary wrap
label_wrap <- function(x) str_wrap(gsub("_", " ", x), width = 30)

# ── Palette (Tableau 10 — colorblind-distinguishable, publication-quality) ─────
ont_colors <- c(
  "BP" = "#4E79A7",   # slate blue
  "MF" = "#59A14F",   # sage green
  "CC" = "#E15759"    # coral red
)

# ── Output naming ──────────────────────────────────────────────────────────────
base_title <- basename(input_file)
base_title <- sub("_TopGo_results_ALL\\.tab$", "", base_title)
base_title <- sub("\\.txt$",                   "", base_title)
out_dir    <- dirname(input_file)    # same folder as input

# ── Shared theme ───────────────────────────────────────────────────────────────
go_theme <- theme_minimal(base_size = 11) +
  theme(
    strip.text          = element_text(face = "bold", size = 11),
    strip.background    = element_rect(fill = "grey95", color = NA),
    axis.text.y         = element_text(size = 8.5, lineheight = 0.9, hjust = 1),
    axis.text.x         = element_text(size = 9),
    axis.title          = element_text(size = 10),
    plot.title          = element_text(face = "bold", size = 13, hjust = 0.5),
    plot.subtitle       = element_text(size = 8.5, hjust = 0.5, color = "grey50"),
    plot.margin         = margin(10, 15, 10, 5),
    legend.title        = element_text(size = 9, face = "bold"),
    legend.text         = element_text(size = 8),
    panel.grid.major.y  = element_blank(),
    panel.grid.minor    = element_blank(),
    panel.border        = element_rect(color = "grey80", fill = NA, linewidth = 0.4)
  )

subtitle <- paste0(
  "p-value correction: ", correction,
  "  \u2022  threshold: p\u202f<\u202f", min_p,
  "  \u2022  fold enrichment\u202f\u2265\u202f", enrich_cutoff
)

# ── Helper: save PDF ──────────────────────────────────────────────────────────
save_plot <- function(p, stem, w, h) {
  pdf_path <- paste0(stem, ".pdf")
  ggsave(pdf_path, plot = p, width = w, height = h, device = "pdf")
  message("  Saved: ", pdf_path)
}

# Estimate wrapped-label height: str_wrap at width 45 can create ~2-line labels,
# so use a taller row allowance than for single-line labels.
n_terms <- nrow(dat)
plot_h  <- max(5, 0.55 * n_terms + 2.5)
plot_w  <- 10

# ══════════════════════════════════════════════════════════════════════════════
# Plot 1 — Horizontal bar chart
# ══════════════════════════════════════════════════════════════════════════════
p_bar <- ggplot(dat, aes(x = logP, y = Term_label, fill = ontology)) +
  geom_col(width = 0.72) +
  geom_vline(
    xintercept = -log10(min_p),
    linetype   = "dashed",
    color      = "grey45",
    linewidth  = 0.4
  ) +
  scale_fill_manual(values = ont_colors, guide = "none") +
  scale_x_continuous(expand = expansion(mult = c(0, 0.08))) +
  scale_y_discrete(labels = label_wrap) +
  facet_grid(ontology ~ ., scales = "free_y", space = "free_y") +
  labs(
    title    = base_title,
    subtitle = subtitle,
    x        = expression(-log[10](italic(p)-value)),
    y        = NULL
  ) +
  go_theme

save_plot(
  p_bar,
  file.path(out_dir, paste0("TopGO_barplot_", base_title)),
  plot_w, plot_h
)

# ══════════════════════════════════════════════════════════════════════════════
# Plot 2 — Dot plot  (fold enrichment × term, size = significant genes,
#                     colour gradient = -log10 p-value)
# ══════════════════════════════════════════════════════════════════════════════
p_dot <- ggplot(dat, aes(
    x     = FoldChange,
    y     = Term_label,
    size  = Significant,
    color = logP
  )) +
  geom_point(alpha = 0.88) +
  scale_color_viridis_c(
    option    = "plasma",
    name      = expression(-log[10](italic(p))),
    direction = 1,
    end       = 0.93
  ) +
  scale_size_continuous(
    name   = "Significant\ngenes",
    range  = c(2.5, 9),
    breaks = function(x) pretty(x, n = 4)
  ) +
  scale_x_continuous(expand = expansion(mult = c(0.05, 0.12))) +
  scale_y_discrete(labels = label_wrap) +
  facet_grid(ontology ~ ., scales = "free_y", space = "free_y") +
  labs(
    title    = base_title,
    subtitle = subtitle,
    x        = "Fold enrichment",
    y        = NULL
  ) +
  go_theme

save_plot(
  p_dot,
  file.path(out_dir, paste0("TopGO_dotplot_", base_title)),
  plot_w, plot_h
)

message("Done: ", base_title)

if (open_plot == 1)
  system(paste("open", shQuote(
    file.path(out_dir, paste0("TopGO_barplot_", base_title, ".pdf"))
  )))
