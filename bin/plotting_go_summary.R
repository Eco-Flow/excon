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

# ── Load & reshape ─────────────────────────────────────────────────────────────
read_go <- function(file, direction) {
  read_tsv(file, show_col_types = FALSE) %>%
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

plot_dat <- dat %>%
  semi_join(sig_terms, by = c("GO_ID", "GO_term", "direction")) %>%
  filter(!is.na(pvalue)) %>%
  mutate(
    neg_log10_p = -log10(pvalue),
    sig         = pvalue < PVAL_CUTOFF,
    GO_label    = str_trunc(GO_term, 45)
  )

term_order <- sig_terms %>%
  group_by(GO_label = str_trunc(GO_term, 45), direction) %>%
  summarise(n_sig = max(n_sig), .groups = "drop") %>%
  arrange(direction, desc(n_sig))

plot_dat <- plot_dat %>%
  mutate(GO_label = factor(GO_label, levels = unique(term_order$GO_label)))

# ── Summary table ──────────────────────────────────────────────────────────────
summary_tbl <- sig_terms %>%
  arrange(direction, desc(n_sig))

cat("\n=== GO Enrichment Summary (p <", PVAL_CUTOFF, ", min", MIN_SPECIES, "species) ===\n")
print(summary_tbl, n = Inf)

write_tsv(summary_tbl, "go_enrichment_summary.tsv")
cat("\nSummary written to go_enrichment_summary.tsv\n")

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

ggsave("go_enrichment_heatmap.pdf", p, width = 16, height = 8, limitsize = FALSE)
ggsave("go_enrichment_heatmap.png", p, width = 16, height = 8, dpi = 150, limitsize = FALSE)
cat("Heatmap saved to go_enrichment_heatmap.pdf / .png\n")

# ── Dot plot ───────────────────────────────────────────────────────────────────
dot_dat <- sig_terms %>%
  mutate(GO_label = factor(str_trunc(GO_term, 45), levels = rev(unique(term_order$GO_label))))

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

ggsave("go_enrichment_dotplot.pdf", p2, width = 10, height = 6, limitsize = FALSE)
ggsave("go_enrichment_dotplot.png", p2, width = 10, height = 6, dpi = 150, limitsize = FALSE)
cat("Dot plot saved to go_enrichment_dotplot.pdf / .png\n")

# ── Aligned figures ────────────────────────────────────────────────────────────
all_sig_ids <- sig_terms %>% distinct(GO_ID, GO_term)

exp_ids  <- sig_terms %>% filter(direction == "expanding")   %>% distinct(GO_ID)
cont_ids <- sig_terms %>% filter(direction == "contracting") %>% distinct(GO_ID)

term_class <- all_sig_ids %>%
  mutate(term_group = case_when(
    GO_ID %in% exp_ids$GO_ID & GO_ID %in% cont_ids$GO_ID ~ "shared",
    GO_ID %in% exp_ids$GO_ID                              ~ "expanding only",
    TRUE                                                   ~ "contracting only"
  ))

label_map <- all_sig_ids %>%
  mutate(GO_label_raw = str_trunc(GO_term, 45)) %>%
  group_by(GO_label_raw) %>%
  mutate(
    GO_label = if (n() > 1) paste0(str_trunc(GO_term, 35), " [", GO_ID, "]") else GO_label_raw
  ) %>%
  ungroup() %>%
  select(GO_ID, GO_term, GO_label)

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
ggsave("go_enrichment_heatmap_aligned.pdf", p3, width = 18, height = h_all, limitsize = FALSE)
ggsave("go_enrichment_heatmap_aligned.png", p3, width = 18, height = h_all, dpi = 150, limitsize = FALSE)
cat("Aligned heatmap (all) saved\n")

# Shared terms only
shared_ids <- term_class %>% filter(term_group == "shared") %>% distinct(GO_ID, GO_term)
if (nrow(shared_ids) > 0) {
  p3b <- make_aligned_heatmap(
    filter(aligned_plot_dat, GO_ID %in% shared_ids$GO_ID),
    " (shared terms only)"
  )
  h_shared <- max(4, nrow(shared_ids) * 0.25 + 2)
  ggsave("go_enrichment_heatmap_aligned_shared.pdf", p3b, width = 18, height = h_shared, limitsize = FALSE)
  ggsave("go_enrichment_heatmap_aligned_shared.png", p3b, width = 18, height = h_shared, dpi = 150, limitsize = FALSE)
  cat("Aligned heatmap (shared) saved\n")
}

# Expanding-only terms
exp_only_ids <- term_class %>% filter(term_group == "expanding only") %>% distinct(GO_ID, GO_term)
if (nrow(exp_only_ids) > 0) {
  p_exp_only <- filter(aligned_plot_dat, GO_ID %in% exp_only_ids$GO_ID, direction == "expanding") %>%
    ggplot(aes(x = species, y = GO_label, fill = neg_log10_p)) +
    geom_tile(colour = "white", linewidth = 0.3) +
    geom_point(data = ~ filter(.x, sig), aes(x = species, y = GO_label),
               shape = 8, size = 1, colour = "black") +
    scale_fill_gradient(low = "lightyellow", high = "#E64B35", na.value = "grey90",
                        name = expression(-log[10](p))) +
    labs(title    = "GO enrichment: expanding-only terms",
         subtitle = paste0("Terms significant in expanding but not contracting; * p < ", PVAL_CUTOFF),
         x = NULL, y = NULL) +
    theme_minimal(base_size = 9) +
    theme(axis.text.x = element_text(angle = 45, hjust = 1, size = 7),
          axis.text.y = element_text(size = 7), panel.grid = element_blank())
  h_exp <- max(4, nrow(exp_only_ids) * 0.25 + 2)
  ggsave("go_enrichment_heatmap_expanding_only.pdf", p_exp_only, width = 14, height = h_exp, limitsize = FALSE)
  ggsave("go_enrichment_heatmap_expanding_only.png", p_exp_only, width = 14, height = h_exp, dpi = 150, limitsize = FALSE)
  cat("Expanding-only heatmap saved\n")
}

# Contracting-only terms
cont_only_ids <- term_class %>% filter(term_group == "contracting only") %>% distinct(GO_ID, GO_term)
if (nrow(cont_only_ids) > 0) {
  p_cont_only <- filter(aligned_plot_dat, GO_ID %in% cont_only_ids$GO_ID, direction == "contracting") %>%
    ggplot(aes(x = species, y = GO_label, fill = neg_log10_p)) +
    geom_tile(colour = "white", linewidth = 0.3) +
    geom_point(data = ~ filter(.x, sig), aes(x = species, y = GO_label),
               shape = 8, size = 1, colour = "black") +
    scale_fill_gradient(low = "lightyellow", high = "#4DBBD5", na.value = "grey90",
                        name = expression(-log[10](p))) +
    labs(title    = "GO enrichment: contracting-only terms",
         subtitle = paste0("Terms significant in contracting but not expanding; * p < ", PVAL_CUTOFF),
         x = NULL, y = NULL) +
    theme_minimal(base_size = 9) +
    theme(axis.text.x = element_text(angle = 45, hjust = 1, size = 7),
          axis.text.y = element_text(size = 7), panel.grid = element_blank())
  h_cont <- max(4, nrow(cont_only_ids) * 0.25 + 2)
  ggsave("go_enrichment_heatmap_contracting_only.pdf", p_cont_only, width = 14, height = h_cont, limitsize = FALSE)
  ggsave("go_enrichment_heatmap_contracting_only.png", p_cont_only, width = 14, height = h_cont, dpi = 150, limitsize = FALSE)
  cat("Contracting-only heatmap saved\n")
}

# ── Aligned dot plot ───────────────────────────────────────────────────────────
full_grid2 <- all_sig_ids %>% crossing(tibble(direction = c("expanding", "contracting")))
aligned_dot_dat <- full_grid2 %>%
  left_join(sig_terms, by = c("GO_ID", "GO_term", "direction")) %>%
  replace_na(list(n_sig = 0)) %>%
  left_join(label_map, by = c("GO_ID", "GO_term")) %>%
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
ggsave("go_enrichment_dotplot_aligned.pdf", p4, width = 12, height = h_dot, limitsize = FALSE)
ggsave("go_enrichment_dotplot_aligned.png", p4, width = 12, height = h_dot, dpi = 150, limitsize = FALSE)
cat("Aligned dot plot saved\n")

cat("\n=== Term overlap summary ===\n")
term_class %>% count(term_group) %>% print()
