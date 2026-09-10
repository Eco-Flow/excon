#!/usr/bin/env Rscript

suppressPackageStartupMessages({
  library(tidyverse)
  library(ggplot2)
})

# ── Parameters ────────────────────────────────────────────────────────────────
PVAL_CUTOFF <- 0.00001
MIN_SPECIES <- 2
POS_FILE    <- "Go_summary_pos.tsv"
NEG_FILE    <- "Go_summary_neg.tsv"

# Write an SVG without hard-depending on svglite (absent from some containers).
# Prefer svglite; else fall back to the cairo-based grDevices::svg(); else skip.
save_svg <- function(plot, file, width, height) {
  if (requireNamespace("svglite", quietly = TRUE)) {
    ggsave(file, plot, width = width, height = height,
           device = svglite::svglite, limitsize = FALSE)
  } else if (capabilities("cairo")) {
    grDevices::svg(file, width = width, height = height)
    print(plot); dev.off()
  } else {
    message("svglite and cairo both unavailable — SVG skipped: ", file)
  }
}

# ── Helper: save PDF + SVG + PNG ──────────────────────────────────────────────
# `stem` may include a subfolder (e.g. "heatmaps/go_enrichment_heatmap"); the
# folder is created on first use so callers don't need a separate mkdir step.
save_plot <- function(p, stem, width, height) {
  dir.create(dirname(stem), recursive = TRUE, showWarnings = FALSE)
  ggsave(paste0(stem, ".pdf"), p, width = width, height = height,
         device = cairo_pdf, limitsize = FALSE)
  save_svg(p, paste0(stem, ".svg"), width, height)
  ggsave(paste0(stem, ".png"), p, width = width, height = height,
         dpi = 150, limitsize = FALSE)
}

# ── Helper: stack two ggplots top/bottom into one PDF + SVG + PNG ─────────────
# Dependency-free (base `grid` only, no patchwork/cowplot), so it works in any
# container that already has ggplot2 without needing an extra package pinned.
# height_top/height_bottom set the relative split, not absolute inches, e.g.
# passing each panel's own natural height keeps the split proportional to how
# many rows it actually has.
save_stacked_plot <- function(p_top, p_bottom, stem, width, height_top, height_bottom) {
  dir.create(dirname(stem), recursive = TRUE, showWarnings = FALSE)
  total_height <- height_top + height_bottom

  draw <- function() {
    grid::grid.newpage()
    grid::pushViewport(grid::viewport(layout = grid::grid.layout(
      nrow = 2, heights = grid::unit(c(height_top, height_bottom), "null")
    )))
    print(p_top,    vp = grid::viewport(layout.pos.row = 1))
    print(p_bottom, vp = grid::viewport(layout.pos.row = 2))
    grid::upViewport(0)
  }

  grDevices::cairo_pdf(paste0(stem, ".pdf"), width = width, height = total_height)
  draw(); grDevices::dev.off()

  if (requireNamespace("svglite", quietly = TRUE)) {
    svglite::svglite(paste0(stem, ".svg"), width = width, height = total_height)
    draw(); grDevices::dev.off()
  } else if (capabilities("cairo")) {
    grDevices::svg(paste0(stem, ".svg"), width = width, height = total_height)
    draw(); grDevices::dev.off()
  } else {
    message("svglite and cairo both unavailable — SVG skipped: ", stem, ".svg")
  }

  grDevices::png(paste0(stem, ".png"), width = width, height = total_height,
                 units = "in", res = 150)
  draw(); grDevices::dev.off()
}

# ── Helper: build a unique, readable label per GO term ────────────────────────
# GO term names are shown close to full length rather than aggressively cut, so
# two different terms don't silently collide into the same-looking label (which
# previously could also merge their data into one heatmap row). Only when two
# terms still tie at `width` characters is the GO ID appended to disambiguate;
# `disambig_width` is kept generous so that fallback still shows real
# distinguishing text, not just an opaque accession number.
make_go_labels <- function(ids_df, width = 70, disambig_width = 55) {
  ids_df %>%
    mutate(GO_label_raw = str_trunc(GO_term, width)) %>%
    group_by(GO_label_raw) %>%
    mutate(
      GO_label = if (n() > 1) paste0(str_trunc(GO_term, disambig_width), " [", GO_ID, "]")
                 else GO_label_raw
    ) %>%
    ungroup() %>%
    select(GO_ID, GO_term, GO_label)
}

# ── Load & reshape ─────────────────────────────────────────────────────────────
read_go <- function(file, direction) {
  read_tsv(file, show_col_types = FALSE, quote = "") %>%
    pivot_longer(
      cols      = -c(GO_ID, GO_term, Count_significant),
      names_to  = "species",
      values_to = "pvalue"
    ) %>%
    mutate(direction = direction)
}

dat <- bind_rows(
  read_go(POS_FILE, "expanding"),
  read_go(NEG_FILE, "contracting")
)

# ── Filter ─────────────────────────────────────────────────────────────────────
sig_terms <- dat %>%
  filter(!is.na(pvalue), pvalue < PVAL_CUTOFF) %>%
  group_by(GO_ID, GO_term, direction) %>%
  summarise(n_sig = n_distinct(species), .groups = "drop") %>%
  filter(n_sig >= MIN_SPECIES)

if (nrow(sig_terms) == 0) {
  cat("\nNo GO terms meet the criteria (p <", PVAL_CUTOFF, ", min", MIN_SPECIES,
      "species) — skipping summary plots.\n")
  quit(status = 0)
}

# Built once, up front, so every plot below uses the same non-colliding labels.
all_sig_ids <- sig_terms %>% distinct(GO_ID, GO_term)
label_map   <- make_go_labels(all_sig_ids)
sig_terms   <- sig_terms %>% left_join(label_map, by = c("GO_ID", "GO_term"))

plot_dat <- dat %>%
  semi_join(sig_terms, by = c("GO_ID", "GO_term", "direction")) %>%
  filter(!is.na(pvalue)) %>%
  left_join(label_map, by = c("GO_ID", "GO_term")) %>%
  mutate(
    neg_log10_p = -log10(pvalue),
    sig         = pvalue < PVAL_CUTOFF
  )

term_order <- sig_terms %>%
  group_by(GO_label, direction) %>%
  summarise(n_sig = max(n_sig), .groups = "drop") %>%
  arrange(direction, desc(n_sig))

plot_dat <- plot_dat %>%
  mutate(GO_label = factor(GO_label, levels = unique(term_order$GO_label)))

# ── Summary table ──────────────────────────────────────────────────────────────
summary_tbl <- sig_terms %>%
  arrange(direction, desc(n_sig))

cat("\n=== GO Enrichment Summary (p <", PVAL_CUTOFF, ", min", MIN_SPECIES, "species) ===\n")
print(summary_tbl, n = Inf)

dir.create("tables", recursive = TRUE, showWarnings = FALSE)
write_tsv(summary_tbl, "tables/go_enrichment_summary.tsv")
cat("\nSummary written to tables/go_enrichment_summary.tsv\n")

# ── Heatmap plot ───────────────────────────────────────────────────────────────
p <- ggplot(plot_dat, aes(x = species, y = GO_label, fill = neg_log10_p)) +
  geom_tile(colour = "white", linewidth = 0.3) +
  geom_point(
    data   = filter(plot_dat, sig),
    aes(x = species, y = GO_label),
    shape = 8, size = 1, colour = "black"
  ) +
  scale_fill_gradient(
    low      = "lightyellow",
    high     = "darkred",
    na.value = "grey90",
    name     = expression(-log[10](p))
  ) +
  facet_wrap(~direction, scales = "free_y") +
  labs(
    title    = "GO enrichment: expanding vs contracting gene families",
    subtitle = paste0("Terms significant (p < ", PVAL_CUTOFF, ") in >= ", MIN_SPECIES,
                      " species shown; * marks p < ", PVAL_CUTOFF),
    x = NULL, y = NULL
  ) +
  theme_minimal(base_size = 9) +
  theme(
    axis.text.x     = element_text(angle = 45, hjust = 1, size = 7),
    axis.text.y     = element_text(size = 7),
    strip.text      = element_text(face = "bold"),
    legend.position = "right",
    panel.grid      = element_blank()
  )

save_plot(p, "heatmaps/go_enrichment_heatmap", 16, 8)
cat("Heatmap saved\n")

# ── Dot plot ───────────────────────────────────────────────────────────────────
dot_dat <- sig_terms %>%
  mutate(GO_label = factor(GO_label, levels = rev(unique(term_order$GO_label))))

p2 <- ggplot(dot_dat, aes(x = n_sig, y = GO_label, colour = direction, size = n_sig)) +
  geom_point(alpha = 0.8) +
  scale_colour_manual(values = c(expanding = "#E64B35", contracting = "#4DBBD5")) +
  scale_size_continuous(range = c(2, 8), guide = "none") +
  labs(
    title  = "Non-redundant GO terms by direction",
    x      = "Number of significant species",
    y      = NULL,
    colour = "Direction"
  ) +
  theme_minimal(base_size = 9) +
  theme(axis.text.y = element_text(size = 7))

save_plot(p2, "dotplots/go_enrichment_dotplot", 10, 6)
cat("Dot plot saved\n")

# ── Aligned figures ────────────────────────────────────────────────────────────
# all_sig_ids and label_map were already built above, right after sig_terms,
# so every plot (including this one) uses the same non-colliding labels.
exp_ids  <- sig_terms %>% filter(direction == "expanding")   %>% distinct(GO_ID)
cont_ids <- sig_terms %>% filter(direction == "contracting") %>% distinct(GO_ID)

term_class <- all_sig_ids %>%
  mutate(term_group = case_when(
    GO_ID %in% exp_ids$GO_ID & GO_ID %in% cont_ids$GO_ID ~ "shared",
    GO_ID %in% exp_ids$GO_ID                              ~ "expanding only",
    TRUE                                                   ~ "contracting only"
  ))

all_species <- dat %>% distinct(species)
complete_grid <- all_sig_ids %>%
  crossing(tibble(direction = c("expanding", "contracting"))) %>%
  crossing(all_species)

shared_order <- sig_terms %>%
  group_by(GO_ID, GO_term) %>%
  summarise(n_total = sum(n_sig), .groups = "drop") %>%
  arrange(desc(n_total)) %>%
  left_join(label_map, by = c("GO_ID", "GO_term")) %>%
  pull(GO_label) %>%
  rev()

aligned_plot_dat <- complete_grid %>%
  left_join(
    dat %>% select(GO_ID, GO_term, species, direction, pvalue),
    by = c("GO_ID", "GO_term", "species", "direction")
  ) %>%
  left_join(label_map, by = c("GO_ID", "GO_term")) %>%
  mutate(
    neg_log10_p = if_else(!is.na(pvalue), -log10(pvalue), NA_real_),
    sig         = !is.na(pvalue) & pvalue < PVAL_CUTOFF,
    GO_label    = factor(GO_label, levels = shared_order),
    direction   = factor(direction, levels = c("expanding", "contracting"))
  )

make_aligned_heatmap <- function(data, title_suffix = "") {
  ggplot(data, aes(x = species, y = GO_label, fill = neg_log10_p)) +
    geom_tile(colour = "white", linewidth = 0.3) +
    geom_point(
      data  = filter(data, sig),
      aes(x = species, y = GO_label),
      shape = 8, size = 1, colour = "black"
    ) +
    scale_fill_gradient(
      low      = "lightyellow",
      high     = "darkred",
      na.value = "grey90",
      name     = expression(-log[10](p))
    ) +
    facet_wrap(~direction, scales = "free_x") +
    labs(
      title    = paste0("GO enrichment: expanding vs contracting", title_suffix),
      subtitle = paste0("Grey = no data in that direction; * p < ", PVAL_CUTOFF),
      x = NULL, y = NULL
    ) +
    theme_minimal(base_size = 9) +
    theme(
      axis.text.x     = element_text(angle = 45, hjust = 1, size = 7),
      axis.text.y     = element_text(size = 7),
      strip.text      = element_text(face = "bold"),
      legend.position = "right",
      panel.grid      = element_blank(),
      panel.spacing   = unit(2, "cm")
    )
}

# All terms
p3 <- make_aligned_heatmap(aligned_plot_dat, " (all terms)")
h_all <- max(4, nrow(all_sig_ids) * 0.25 + 2)
save_plot(p3, "heatmaps/go_enrichment_heatmap_aligned", 18, h_all)
cat("Aligned heatmap (all) saved\n")

# Shared terms only
shared_ids <- term_class %>% filter(term_group == "shared") %>% distinct(GO_ID, GO_term)
if (nrow(shared_ids) > 0) {
  p3b <- make_aligned_heatmap(
    filter(aligned_plot_dat, GO_ID %in% shared_ids$GO_ID),
    " (shared terms only)"
  )
  h_shared <- max(4, nrow(shared_ids) * 0.25 + 2)
  save_plot(p3b, "heatmaps/go_enrichment_heatmap_aligned_shared", 18, h_shared)
  cat("Aligned heatmap (shared) saved\n")
}

# Expanding-only terms
p_exp_only <- NULL
h_exp      <- NULL
exp_only_ids <- term_class %>% filter(term_group == "expanding only") %>% distinct(GO_ID, GO_term)
if (nrow(exp_only_ids) > 0) {
  p_exp_only <- filter(aligned_plot_dat, GO_ID %in% exp_only_ids$GO_ID,
                       direction == "expanding") %>%
    ggplot(aes(x = species, y = GO_label, fill = neg_log10_p)) +
    geom_tile(colour = "white", linewidth = 0.3) +
    geom_point(data = ~ filter(.x, sig), aes(x = species, y = GO_label),
               shape = 8, size = 1, colour = "black") +
    scale_fill_gradient(low = "lightyellow", high = "#E64B35", na.value = "grey90",
                        name = expression(-log[10](p))) +
    labs(title    = "GO enrichment: expanding-only terms",
         subtitle = paste0("Terms significant in expanding but not contracting; * p < ",
                           PVAL_CUTOFF),
         x = NULL, y = NULL) +
    theme_minimal(base_size = 9) +
    theme(axis.text.x = element_text(angle = 45, hjust = 1, size = 7),
          axis.text.y = element_text(size = 7), panel.grid = element_blank())
  h_exp <- max(4, nrow(exp_only_ids) * 0.25 + 2)
  save_plot(p_exp_only, "heatmaps/go_enrichment_heatmap_expanding_only", 14, h_exp)
  cat("Expanding-only heatmap saved\n")
}

# Contracting-only terms
p_cont_only <- NULL
h_cont      <- NULL
cont_only_ids <- term_class %>% filter(term_group == "contracting only") %>% distinct(GO_ID, GO_term)
if (nrow(cont_only_ids) > 0) {
  p_cont_only <- filter(aligned_plot_dat, GO_ID %in% cont_only_ids$GO_ID,
                        direction == "contracting") %>%
    ggplot(aes(x = species, y = GO_label, fill = neg_log10_p)) +
    geom_tile(colour = "white", linewidth = 0.3) +
    geom_point(data = ~ filter(.x, sig), aes(x = species, y = GO_label),
               shape = 8, size = 1, colour = "black") +
    scale_fill_gradient(low = "lightyellow", high = "#4DBBD5", na.value = "grey90",
                        name = expression(-log[10](p))) +
    labs(title    = "GO enrichment: contracting-only terms",
         subtitle = paste0("Terms significant in contracting but not expanding; * p < ",
                           PVAL_CUTOFF),
         x = NULL, y = NULL) +
    theme_minimal(base_size = 9) +
    theme(axis.text.x = element_text(angle = 45, hjust = 1, size = 7),
          axis.text.y = element_text(size = 7), panel.grid = element_blank())
  h_cont <- max(4, nrow(cont_only_ids) * 0.25 + 2)
  save_plot(p_cont_only, "heatmaps/go_enrichment_heatmap_contracting_only", 14, h_cont)
  cat("Contracting-only heatmap saved\n")
}

# Expanding-only + contracting-only stacked into one combined figure (top/bottom),
# matching the layout used for the paper's Figure 3. Each panel keeps its own
# height proportional to its row count, same as the two figures do standalone;
# only built when both panels actually exist.
if (!is.null(p_exp_only) && !is.null(p_cont_only)) {
  save_stacked_plot(
    p_exp_only, p_cont_only,
    "heatmaps/go_enrichment_heatmap_expanding_contracting_stacked",
    width = 14, height_top = h_exp, height_bottom = h_cont
  )
  cat("Expanding+contracting stacked heatmap saved\n")
}

# ── Aligned dot plot ───────────────────────────────────────────────────────────
full_grid2 <- all_sig_ids %>% crossing(tibble(direction = c("expanding", "contracting")))
aligned_dot_dat <- full_grid2 %>%
  left_join(sig_terms, by = c("GO_ID", "GO_term", "direction")) %>%
  replace_na(list(n_sig = 0)) %>%
  left_join(term_class %>% select(GO_ID, term_group), by = "GO_ID") %>%
  mutate(
    GO_label  = factor(GO_label, levels = shared_order),
    direction = factor(direction, levels = c("expanding", "contracting")),
    present   = n_sig >= MIN_SPECIES
  )

p4 <- ggplot(aligned_dot_dat, aes(x = n_sig, y = GO_label, colour = direction, size = n_sig)) +
  geom_point(data = filter(aligned_dot_dat, present), alpha = 0.8) +
  scale_colour_manual(values = c(expanding = "#E64B35", contracting = "#4DBBD5")) +
  scale_size_continuous(range = c(2, 8), guide = "none") +
  facet_wrap(~direction) +
  labs(
    title  = "Non-redundant GO terms by direction (aligned)",
    x      = "Number of significant species",
    y      = NULL,
    colour = "Direction"
  ) +
  theme_minimal(base_size = 9) +
  theme(
    axis.text.y   = element_text(size = 7),
    strip.text    = element_text(face = "bold"),
    panel.spacing = unit(2, "cm")
  )

h_dot <- max(4, nrow(all_sig_ids) * 0.2 + 2)
save_plot(p4, "dotplots/go_enrichment_dotplot_aligned", 12, h_dot)
cat("Aligned dot plot saved\n")

cat("\n=== Term overlap summary ===\n")
term_class %>% count(term_group) %>% print()

# ── README ─────────────────────────────────────────────────────────────────────
readme_lines <- c(
  "GO enrichment plots — output guide",
  "===================================",
  "",
  paste0("Terms shown: significant at p < ", PVAL_CUTOFF, " in >= ", MIN_SPECIES, " species."),
  "Every figure is written as .pdf (vector, for print), .svg (vector, editable",
  "in Inkscape/Illustrator) and .png (quick preview).",
  "",
  "heatmaps/",
  "  go_enrichment_heatmap                          Original two-panel heatmap",
  "                                                  (expanding | contracting),",
  "                                                  every significant term shown",
  "                                                  once per panel independently.",
  "  go_enrichment_heatmap_aligned                   Same data, but expanding and",
  "                                                  contracting panels share one",
  "                                                  term ordering/axis, so the",
  "                                                  same row is the same GO term",
  "                                                  in both panels. All terms.",
  "  go_enrichment_heatmap_aligned_shared             Aligned version, restricted to",
  "                                                  terms significant in BOTH",
  "                                                  directions.",
  "  go_enrichment_heatmap_expanding_only             Aligned version, restricted to",
  "                                                  terms significant only among",
  "                                                  expanding families.",
  "  go_enrichment_heatmap_contracting_only           Aligned version, restricted to",
  "                                                  terms significant only among",
  "                                                  contracting families.",
  "  go_enrichment_heatmap_expanding_contracting_stacked",
  "                                                  Expanding-only (top) and",
  "                                                  contracting-only (bottom)",
  "                                                  combined into one figure —",
  "                                                  used as Figure 3 in the paper.",
  "",
  "dotplots/",
  "  go_enrichment_dotplot                           Non-redundant term list, one",
  "                                                  dot per term (size/position =",
  "                                                  number of significant species),",
  "                                                  expanding and contracting on",
  "                                                  the same axis.",
  "  go_enrichment_dotplot_aligned                   Same, but faceted expanding vs",
  "                                                  contracting with a shared term",
  "                                                  ordering, matching the aligned",
  "                                                  heatmaps above.",
  "",
  "tables/",
  "  go_enrichment_summary.tsv                       One row per significant",
  "                                                  (GO term, direction) pair:",
  "                                                  GO_ID, GO_term, direction,",
  "                                                  n_sig (number of species where",
  "                                                  this term was significant),",
  "                                                  GO_label (the plot axis label).",
  "",
  "Generated by bin/plotting_go_summary.R."
)
writeLines(readme_lines, "README.md")
cat("\nREADME.md written\n")
