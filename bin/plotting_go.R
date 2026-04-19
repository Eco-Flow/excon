#!/usr/bin/env Rscript

suppressPackageStartupMessages({
  library(dplyr)
  library(tidyr)
  library(ggplot2)
})

p_threshold   <- 0.05
shared_cutoff <- 2   # GO term must be significant in >= this many species

make_heatmap <- function(tsv_file, outfile, title, exclude_nodes = FALSE) {
    if (!file.exists(tsv_file)) {
        message(sprintf("File not found: %s — skipping", tsv_file))
        return(invisible(NULL))
    }

    df <- read.table(tsv_file, header = TRUE, sep = "\t", check.names = FALSE,
                     stringsAsFactors = FALSE)

    df_shared <- df[df$Count_significant >= shared_cutoff, ]

    if (nrow(df_shared) == 0) {
        pdf(outfile, width = 7, height = 5)
        plot.new()
        text(0.5, 0.5,
             sprintf("No GO terms shared across >= %d species", shared_cutoff),
             cex = 1.2)
        dev.off()
        message(sprintf("No shared terms — blank plot written to %s", outfile))
        return(invisible(NULL))
    }

    df_shared$Term_Label <- paste0(substr(df_shared$GO_term, 1, 45),
                                   " (", df_shared$GO_ID, ")")

    species_cols <- setdiff(colnames(df_shared),
                            c("GO_ID", "GO_term", "Count_significant", "Term_Label"))
    if (exclude_nodes) {
        species_cols <- species_cols[!grepl("Node", species_cols)]
    }
    if (length(species_cols) == 0) {
        message(sprintf("No species columns remain after filtering — skipping %s", outfile))
        return(invisible(NULL))
    }

    term_order <- df_shared %>%
        arrange(desc(Count_significant)) %>%
        pull(Term_Label)

    long_df <- df_shared %>%
        select(Term_Label, all_of(species_cols)) %>%
        pivot_longer(cols      = all_of(species_cols),
                     names_to  = "Species",
                     values_to = "pvalue") %>%
        mutate(
            pvalue      = ifelse(is.na(pvalue), 1, pvalue),
            neg_log10_p = -log10(pvalue),
            sig_label   = ifelse(pvalue < p_threshold, "*", ""),
            Term_Label  = factor(Term_Label, levels = rev(term_order))
        )

    n_terms   <- length(unique(long_df$Term_Label))
    n_species <- length(species_cols)
    plot_h    <- max(5, n_terms * 0.35 + 2)
    plot_w    <- max(7, n_species * 0.9 + 4)

    p <- ggplot(long_df, aes(x = Species, y = Term_Label, fill = neg_log10_p)) +
        geom_tile(color = "white", linewidth = 0.4) +
        geom_text(aes(label = sig_label), size = 4, vjust = 0.75) +
        scale_fill_gradient(low = "white", high = "#D95F0E",
                            name = "-log10(p)") +
        theme_minimal(base_size = 10) +
        theme(
            axis.text.x = element_text(angle = 45, hjust = 1),
            axis.text.y = element_text(size  = 8),
            panel.grid  = element_blank()
        ) +
        labs(
            title    = title,
            subtitle = sprintf("GO terms significant in >= %d species  |  * p < %.2f",
                               shared_cutoff, p_threshold),
            x = NULL, y = NULL
        )

    ggsave(outfile, p, width = plot_w, height = plot_h)
    message(sprintf("Saved: %s  (%d terms x %d species)", outfile, n_terms, n_species))
}

make_heatmap("Go_summary_pos.tsv",
             "Go_summary_pos.pdf",
             "Gene Family Expansions — Shared GO Terms")

make_heatmap("Go_summary_pos.tsv",
             "Go_summary_pos_noNode.pdf",
             "Gene Family Expansions — Shared GO Terms (species only)",
             exclude_nodes = TRUE)

make_heatmap("Go_summary_neg.tsv",
             "Go_summary_neg.pdf",
             "Gene Family Contractions — Shared GO Terms")

make_heatmap("Go_summary_neg.tsv",
             "Go_summary_neg_noNode.pdf",
             "Gene Family Contractions — Shared GO Terms (species only)",
             exclude_nodes = TRUE)
