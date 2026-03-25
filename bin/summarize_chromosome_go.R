#!/usr/bin/Rscript

# ============================================================================
# GO ENRICHMENT SUMMARY (CLEAN VERSION)
# ============================================================================

suppressPackageStartupMessages({
  library(dplyr)
  library(tidyr)
  library(ggplot2)
  library(optparse)
})

# ============================================================================
# ARGUMENT PARSING
# ============================================================================

option_list <- list(
  make_option(c("--bonferroni"), type="double", default=0.01,
              help="Bonferroni threshold [default %default]"),
  make_option(c("--min_annotated"), type="integer", default=10,
              help="Minimum annotated genes [default %default]"),
  make_option(c("--top_n_hits"), type="integer", default=5,
              help="Top hits per chromosome [default %default]")
)

opt <- parse_args(OptionParser(option_list=option_list))

bonferroni_threshold <- opt$bonferroni
min_annotated <- opt$min_annotated
top_n_hits <- opt$top_n_hits

# ============================================================================
# 2. LOAD DATA
# ============================================================================

go_files <- list.files(path = go_results_dir,
                       pattern = file_pattern,
                       full.names = TRUE)

cat(sprintf("Found %d GO files\n", length(go_files)))

if (length(go_files) == 0) {
  stop("No GO files found")
}

all_go_results <- lapply(go_files, function(f) {
  chrom <- gsub(".*_([0-9]+)_res\\.tab", "\\1", basename(f))

  df <- read.table(f, header = TRUE, sep = "\t",
                   quote = "", comment.char = "",
                   stringsAsFactors = FALSE)

  df$Chromosome <- paste0("Scaffold_", chrom)
  df$Chromosome_Num <- as.numeric(chrom)

  df
}) %>% bind_rows()

cat(sprintf("Total GO terms: %d\n", nrow(all_go_results)))

# ============================================================================
# 3. FILTER (KEY FIX HERE)
# ============================================================================

significant_go <- all_go_results %>%
  filter(
    bonferroni < bonferroni_threshold,
    Annotated >= min_annotated
  ) %>%
  arrange(Chromosome_Num, bonferroni)

cat(sprintf("Significant terms (filtered): %d\n", nrow(significant_go)))

# ============================================================================
# 4. SUMMARY PER CHROMOSOME
# ============================================================================

sig_per_chrom <- significant_go %>%
  count(Chromosome, name = "N_Significant") %>%
  arrange(desc(N_Significant))

print(sig_per_chrom)

# ============================================================================
# 5. TOP HITS
# ============================================================================

top_hits <- significant_go %>%
  group_by(Chromosome) %>%
  slice_min(order_by = bonferroni, n = top_n_hits) %>%
  ungroup()

# ============================================================================
# 6. PLOT 1 — SIGNIFICANT TERMS PER CHROMOSOME
# ============================================================================

p1 <- ggplot(sig_per_chrom,
             aes(x = reorder(Chromosome, N_Significant),
                 y = N_Significant)) +
  geom_col(fill = "#377EB8") +
  coord_flip() +
  theme_minimal() +
  labs(
    title = "Significant GO Terms per Chromosome",
    subtitle = sprintf("Bonferroni < %.2f | Annotated >= %d",
                       bonferroni_threshold, min_annotated),
    x = "Chromosome",
    y = "Count"
  )

ggsave("go_summary_per_chromosome.pdf", p1, width = 8, height = 6)

# ============================================================================
# 7. HEATMAP DATA
# ============================================================================

heatmap_data <- top_hits %>%
  mutate(
    Term_Label = paste0(substr(Term, 1, 40), " (", GO.ID, ")"),
    Neg_log10_P = -log10(bonferroni),
    Percent_Significant = 100 * Significant / Annotated
  )

# ============================================================================
# 8. HEATMAP
# ============================================================================

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
    title = "GO Enrichment Heatmap",
    subtitle = sprintf("Bonferroni < %.2e | Annotated >= %d",
                       bonferroni_threshold, min_annotated),
    x = "Chromosome",
    y = "GO Term"
  )

ggsave("go_heatmap.pdf", p2,
       width = 12,
       height = max(6, nrow(heatmap_data) * 0.3))

# ============================================================================
# 9. ONTOLOGY SUMMARY
# ============================================================================

ontology_summary <- significant_go %>%
  count(Chromosome, ontology)

p3 <- ggplot(ontology_summary,
             aes(x = Chromosome,
                 y = n,
                 fill = ontology)) +
  geom_col() +
  theme_minimal() +
  theme(axis.text.x = element_text(angle = 45, hjust = 1)) +
  labs(
    title = "GO Ontology Distribution",
    subtitle = sprintf("Bonferroni < %.2f | Annotated >= %d",
                       bonferroni_threshold, min_annotated),
    x = "Chromosome",
    y = "Count"
  )

ggsave("go_summary_ontology.pdf", p3, width = 10, height = 6)

# ============================================================================
# 10. EXPORT TABLES
# ============================================================================

# Full filtered set
write.csv(significant_go,
          "go_significant_filtered.csv",
          row.names = FALSE)

# Top hits
write.csv(top_hits,
          "go_top_hits.csv",
          row.names = FALSE)

# Summary stats
summary_stats <- all_go_results %>%
  group_by(Chromosome) %>%
  summarise(
    Total = n(),
    Significant = sum(bonferroni < bonferroni_threshold &
                      Annotated >= min_annotated),
    .groups = "drop"
  )

write.csv(summary_stats,
          "go_summary_stats.csv",
          row.names = FALSE)

# ============================================================================
# DONE
# ============================================================================

cat("\n=== DONE ===\n")
cat("Applied filters:\n")
cat(sprintf("  Bonferroni < %.2f\n", bonferroni_threshold))
cat(sprintf("  Min annotated >= %d\n", min_annotated))
cat(sprintf("  Top hits per chromosome = %d\n\n", top_n_hits))
