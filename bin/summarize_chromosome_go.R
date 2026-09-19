#!/usr/bin/env Rscript

# ============================================================================
# GO ENRICHMENT SUMMARY
# ============================================================================

suppressPackageStartupMessages({
  library(dplyr)
  library(tidyr)
  library(ggplot2)
})

# ============================================================================
# ARGUMENT PARSING (base R — no optparse dependency)
# ============================================================================

args <- commandArgs(trailingOnly = TRUE)

parse_arg <- function(args, flag, default = NULL) {
  idx <- which(args == flag)
  if (length(idx) && idx < length(args)) args[idx + 1] else default
}

go_results_dir       <- parse_arg(args, "--input")
bonferroni_threshold <- as.double(parse_arg(args, "--bonferroni", 0.01))
min_annotated        <- as.integer(parse_arg(args, "--min_annotated", 10))
top_n_hits           <- as.integer(parse_arg(args, "--top_n_hits", 5))
max_chroms           <- as.integer(parse_arg(args, "--max_chroms", 15))

if (is.null(go_results_dir)) stop("--input is required")
if (!dir.exists(go_results_dir)) stop(sprintf("Input directory not found: %s", go_results_dir))

# Prefix all outputs with the input directory name so filtered/unfiltered
# and per-species outputs never collide when collected downstream
prefix <- basename(go_results_dir)

cat(sprintf("Input directory : %s\n", go_results_dir))
cat(sprintf("Output prefix   : %s\n", prefix))

# ============================================================================
# LOAD DATA
# ============================================================================

go_files <- list.files(path = go_results_dir,
                       pattern = "_res\\.tab$",
                       full.names = TRUE)

cat(sprintf("Found %d GO result files\n", length(go_files)))

if (length(go_files) == 0) {
  stop(sprintf("No *_res.tab files found in: %s", go_results_dir))
}

all_go_results <- lapply(go_files, function(f) {
  chrom <- gsub(".*_([0-9]+)_res\\.tab$", "\\1", basename(f))

  df <- read.table(f, header = TRUE, sep = "\t",
                   quote = "", comment.char = "",
                   stringsAsFactors = FALSE)

  df$Chromosome     <- paste0("Scaffold_", chrom)
  df$Chromosome_Num <- as.numeric(chrom)
  df
}) %>% bind_rows()

cat(sprintf("Total GO terms loaded: %d\n", nrow(all_go_results)))

# ============================================================================
# FILTER
# ============================================================================

significant_go <- all_go_results %>%
  filter(
    bonferroni < bonferroni_threshold,
    Annotated  >= min_annotated
  ) %>%
  arrange(Chromosome_Num, bonferroni)

cat(sprintf("Significant terms after filtering: %d\n", nrow(significant_go)))

# ============================================================================
# SUMMARY PER CHROMOSOME
# ============================================================================

sig_per_chrom <- significant_go %>%
  count(Chromosome, name = "N_Significant") %>%
  arrange(desc(N_Significant))

print(sig_per_chrom)

# Plots below are capped to the --max_chroms most significant chromosomes.
# Without this, a genome with many scaffolds (a fragmented assembly, or just a
# species with a large karyotype — birds routinely have 30-80 chromosomes)
# produces a heatmap tall enough to exceed ggsave()'s 50-inch sanity limit and
# fail outright; even short of that ceiling, a plot with dozens of chromosomes
# stops being readable well before it stops being drawable. This only trims
# the *plots* — every chromosome's full results are still in the CSV outputs.
n_chroms_total <- nrow(sig_per_chrom)
chroms_to_plot <- sig_per_chrom %>%
  slice_head(n = max_chroms) %>%
  pull(Chromosome)
n_chroms_dropped <- n_chroms_total - length(chroms_to_plot)

if (n_chroms_dropped > 0) {
  cat(sprintf(
    "Plots show the top %d of %d chromosomes by significant-term count (%d omitted from plots only — full data still in the CSV outputs; adjust with --max_chroms)\n",
    length(chroms_to_plot), n_chroms_total, n_chroms_dropped
  ))
}

plot_subtitle_suffix <- if (n_chroms_dropped > 0) {
  sprintf(" | top %d of %d chromosomes shown", length(chroms_to_plot), n_chroms_total)
} else {
  ""
}

# ============================================================================
# TOP HITS
# ============================================================================

top_hits <- significant_go %>%
  group_by(Chromosome) %>%
  slice_min(order_by = bonferroni, n = top_n_hits) %>%
  ungroup()

# ============================================================================
# PLOT 1 — SIGNIFICANT TERMS PER CHROMOSOME
# ============================================================================

p1 <- ggplot(sig_per_chrom %>% filter(Chromosome %in% chroms_to_plot),
             aes(x = reorder(Chromosome, N_Significant),
                 y = N_Significant)) +
  geom_col(fill = "#377EB8") +
  coord_flip() +
  theme_minimal() +
  labs(
    title    = "Significant GO Terms per Chromosome",
    subtitle = paste0(sprintf("Bonferroni < %.2f | Annotated >= %d",
                       bonferroni_threshold, min_annotated),
                       plot_subtitle_suffix),
    x = "Chromosome",
    y = "Count"
  )

ggsave(paste0(prefix, "_go_summary_per_chromosome.pdf"), p1,
       width = 8, height = 6)

# ============================================================================
# PLOT 2 — HEATMAP
# ============================================================================

heatmap_data <- top_hits %>%
  filter(Chromosome %in% chroms_to_plot) %>%
  mutate(
    Term_Label           = paste0(substr(Term, 1, 40), " (", GO.ID, ")"),
    Neg_log10_P          = -log10(bonferroni),
    Percent_Significant  = 100 * Significant / Annotated
  )

p2 <- ggplot(heatmap_data,
             aes(x = Chromosome,
                 y = Term_Label,
                 fill = Neg_log10_P)) +
  geom_tile(color = "white") +
  geom_text(aes(label = sprintf("%.1f%%", Percent_Significant)),
            size = 2.5) +
  scale_fill_gradient(low = "#FFF7BC", high = "#D95F0E") +
  theme_minimal(base_size = 10) +
  theme(axis.text.x = element_text(angle = 45, hjust = 1)) +
  labs(
    title    = "GO Enrichment Heatmap",
    subtitle = paste0(sprintf("Bonferroni < %.2e | Annotated >= %d",
                       bonferroni_threshold, min_annotated),
                       plot_subtitle_suffix),
    x = "Chromosome",
    y = "GO Term"
  )

ggsave(paste0(prefix, "_go_heatmap.pdf"), p2,
       width  = 12,
       # Capped at 48in regardless of --max_chroms/--top_n_hits, so a generous
       # setting densifies the plot instead of exceeding ggsave()'s 50in limit.
       height = min(48, max(6, nrow(heatmap_data) * 0.3)))

# ============================================================================
# PLOT 3 — ONTOLOGY DISTRIBUTION
# ============================================================================

ontology_summary <- significant_go %>%
  count(Chromosome, ontology)

p3 <- ggplot(ontology_summary %>% filter(Chromosome %in% chroms_to_plot),
             aes(x = Chromosome,
                 y = n,
                 fill = ontology)) +
  geom_col() +
  theme_minimal() +
  theme(axis.text.x = element_text(angle = 45, hjust = 1)) +
  labs(
    title    = "GO Ontology Distribution",
    subtitle = paste0(sprintf("Bonferroni < %.2f | Annotated >= %d",
                       bonferroni_threshold, min_annotated),
                       plot_subtitle_suffix),
    x = "Chromosome",
    y = "Count"
  )

ggsave(paste0(prefix, "_go_summary_ontology.pdf"), p3,
       width = 10, height = 6)

# ============================================================================
# EXPORT TABLES
# ============================================================================

write.csv(significant_go,
          paste0(prefix, "_go_significant_filtered.csv"),
          row.names = FALSE)

write.csv(top_hits,
          paste0(prefix, "_go_top_hits.csv"),
          row.names = FALSE)

summary_stats <- all_go_results %>%
  group_by(Chromosome) %>%
  summarise(
    Total       = n(),
    Significant = sum(bonferroni < bonferroni_threshold &
                      Annotated  >= min_annotated),
    .groups = "drop"
  )

write.csv(summary_stats,
          paste0(prefix, "_go_summary_stats.csv"),
          row.names = FALSE)

# ============================================================================
# DONE
# ============================================================================

cat("\n=== DONE ===\n")
cat(sprintf("  Output prefix          : %s\n", prefix))
cat(sprintf("  Bonferroni threshold   : %.2f\n", bonferroni_threshold))
cat(sprintf("  Min annotated          : %d\n",   min_annotated))
cat(sprintf("  Top hits per chromosome: %d\n", top_n_hits))
cat(sprintf("  Max chromosomes plotted: %d\n\n", max_chroms))
