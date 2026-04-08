#!/usr/bin/env Rscript
# =============================================================================
# plot_benchmark.R
# Reads benchmark_metrics.tsv + benchmark_per_process.tsv produced by
# collect_metrics.py and generates publication-quality scaling figures.
#
# Usage:
#   Rscript benchmark/plot_benchmark.R \
#       --metrics  benchmark/results/benchmark_metrics.tsv \
#       --perproc  benchmark/results/benchmark_per_process.tsv \
#       --outdir   benchmark/results/figures
# =============================================================================

suppressPackageStartupMessages({
  library(ggplot2)
  library(dplyr)
  library(tidyr)
  library(scales)
  library(patchwork)
})

# ---- Argument parsing -------------------------------------------------------
args <- commandArgs(trailingOnly = TRUE)

get_arg <- function(flag, default = NULL) {
  i <- which(args == flag)
  if (length(i) && length(args) >= i + 1) args[i + 1] else default
}

metrics_file <- get_arg("--metrics", "benchmark/results/benchmark_metrics.tsv")
perproc_file <- get_arg("--perproc", "benchmark/results/benchmark_per_process.tsv")
outdir       <- get_arg("--outdir",  "benchmark/results/figures")

if (!file.exists(metrics_file)) stop("Cannot find --metrics file: ", metrics_file)
if (!file.exists(perproc_file)) stop("Cannot find --perproc file: ", perproc_file)
dir.create(outdir, recursive = TRUE, showWarnings = FALSE)

# ---- Load data --------------------------------------------------------------
metrics <- read.delim(metrics_file, stringsAsFactors = FALSE) |>
  mutate(
    genome_size = factor(genome_size, levels = c("bacteria", "insect", "mammal")),
    quality     = factor(quality,     levels = c("contiguous", "fragmented")),
    phylogeny   = factor(phylogeny,   levels = c("close", "diverse")),
    n_species   = as.integer(n_species)
  )

perproc <- read.delim(perproc_file, stringsAsFactors = FALSE) |>
  mutate(
    genome_size = factor(genome_size, levels = c("bacteria", "insect", "mammal")),
    quality     = factor(quality,     levels = c("contiguous", "fragmented")),
    phylogeny   = factor(phylogeny,   levels = c("close", "diverse")),
    n_species   = as.integer(n_species),
    wall_time_min = wall_time_s / 60
  )

# ---- Shared theme -----------------------------------------------------------
theme_bench <- function() {
  theme_bw(base_size = 11) +
  theme(
    strip.background  = element_rect(fill = "grey92", colour = "grey70"),
    strip.text        = element_text(face = "bold"),
    legend.position   = "bottom",
    panel.grid.minor  = element_blank()
  )
}

genome_colours <- c(bacteria = "#4DAF4A", insect = "#377EB8", mammal = "#E41A1C")
quality_ltypes  <- c(contiguous = "solid", fragmented = "dashed")

# =============================================================================
# Figure 1 — Overall scaling: wall time + peak RAM
# =============================================================================
p_wall <- ggplot(metrics,
    aes(x = n_species, y = total_wall_time_min,
        colour = genome_size, linetype = quality, shape = quality)) +
  geom_point(size = 2.5) +
  geom_line() +
  facet_wrap(~phylogeny, labeller = label_both) +
  scale_colour_manual(values = genome_colours, name = "Genome size") +
  scale_linetype_manual(values = quality_ltypes, name = "Assembly quality") +
  scale_shape_manual(values = c(contiguous = 16, fragmented = 1), name = "Assembly quality") +
  scale_x_continuous(breaks = unique(metrics$n_species)) +
  scale_y_continuous(labels = comma) +
  labs(x = "Number of genomes", y = "Total wall time (min)") +
  theme_bench()

p_ram <- ggplot(metrics,
    aes(x = n_species, y = peak_ram_gb,
        colour = genome_size, linetype = quality, shape = quality)) +
  geom_point(size = 2.5) +
  geom_line() +
  facet_wrap(~phylogeny, labeller = label_both) +
  scale_colour_manual(values = genome_colours, name = "Genome size") +
  scale_linetype_manual(values = quality_ltypes, name = "Assembly quality") +
  scale_shape_manual(values = c(contiguous = 16, fragmented = 1), name = "Assembly quality") +
  scale_x_continuous(breaks = unique(metrics$n_species)) +
  labs(x = "Number of genomes", y = "Peak RAM (GB)") +
  theme_bench()

fig1 <- p_wall / p_ram +
  plot_annotation(
    title = "EXCON pipeline scaling",
    tag_levels = "A"
  ) +
  plot_layout(guides = "collect") &
  theme(legend.position = "bottom")

ggsave(file.path(outdir, "fig1_overall_scaling.pdf"), fig1,
       width = 8, height = 7, useDingbats = FALSE)
ggsave(file.path(outdir, "fig1_overall_scaling.png"), fig1,
       width = 8, height = 7, dpi = 300)
message("Saved: fig1_overall_scaling")

# =============================================================================
# Figure 2 — Per-module wall-time scaling
# =============================================================================

# Key pipeline stages in display order
stage_order <- c(
  "NCBIGENOMEDOWNLOAD", "AGAT_SPKEEPLONGESTISOFORM",
  "ORTHOFINDER_CAFE", "ORTHOFINDER_V2_CAFE",
  "RESCALE_TREE", "CAFE_PREP",
  "CAFE_RUN_K", "CAFE_SELECT_K",
  "CAFE_RUN_BEST", "CAFE_MODEL_COMPARE", "CAFE_PLOT"
)

perproc_key <- perproc |>
  filter(process %in% stage_order) |>
  mutate(process = factor(process, levels = stage_order))

p_modules <- ggplot(perproc_key,
    aes(x = n_species, y = wall_time_min,
        colour = genome_size, linetype = quality)) +
  geom_point(size = 1.8, aes(shape = quality)) +
  geom_line() +
  facet_wrap(~process, scales = "free_y", ncol = 3) +
  scale_colour_manual(values = genome_colours, name = "Genome size") +
  scale_linetype_manual(values = quality_ltypes, name = "Assembly quality") +
  scale_shape_manual(values = c(contiguous = 16, fragmented = 1), name = "Assembly quality") +
  scale_x_continuous(breaks = unique(perproc_key$n_species)) +
  labs(
    x = "Number of genomes",
    y = "Wall time (min)",
    title = "Per-module scaling"
  ) +
  theme_bench() +
  theme(axis.text.x = element_text(size = 8))

ggsave(file.path(outdir, "fig2_per_module_scaling.pdf"), p_modules,
       width = 10, height = 9, useDingbats = FALSE)
ggsave(file.path(outdir, "fig2_per_module_scaling.png"), p_modules,
       width = 10, height = 9, dpi = 300)
message("Saved: fig2_per_module_scaling")

# =============================================================================
# Figure 3 — CPU efficiency + throughput
# =============================================================================
p_eff <- ggplot(metrics,
    aes(x = n_species, y = cpu_efficiency * 100,
        colour = genome_size, linetype = quality, shape = quality)) +
  geom_point(size = 2.5) +
  geom_line() +
  facet_wrap(~phylogeny, labeller = label_both) +
  scale_colour_manual(values = genome_colours, name = "Genome size") +
  scale_linetype_manual(values = quality_ltypes, name = "Assembly quality") +
  scale_shape_manual(values = c(contiguous = 16, fragmented = 1), name = "Assembly quality") +
  scale_x_continuous(breaks = unique(metrics$n_species)) +
  scale_y_continuous(limits = c(0, NA), labels = function(x) paste0(x, "%")) +
  labs(x = "Number of genomes", y = "CPU efficiency") +
  theme_bench()

p_tput <- ggplot(metrics,
    aes(x = n_species, y = genomes_per_hour,
        colour = genome_size, linetype = quality, shape = quality)) +
  geom_point(size = 2.5) +
  geom_line() +
  facet_wrap(~phylogeny, labeller = label_both) +
  scale_colour_manual(values = genome_colours, name = "Genome size") +
  scale_linetype_manual(values = quality_ltypes, name = "Assembly quality") +
  scale_shape_manual(values = c(contiguous = 16, fragmented = 1), name = "Assembly quality") +
  scale_x_continuous(breaks = unique(metrics$n_species)) +
  labs(x = "Number of genomes", y = "Throughput (genomes / hour)") +
  theme_bench()

fig3 <- p_eff / p_tput +
  plot_annotation(
    title = "Pipeline efficiency and throughput",
    tag_levels = "A"
  ) +
  plot_layout(guides = "collect") &
  theme(legend.position = "bottom")

ggsave(file.path(outdir, "fig3_efficiency.pdf"), fig3,
       width = 8, height = 7, useDingbats = FALSE)
ggsave(file.path(outdir, "fig3_efficiency.png"), fig3,
       width = 8, height = 7, dpi = 300)
message("Saved: fig3_efficiency")

# =============================================================================
# Figure 4 — Stacked bar: time composition at maximum dataset size
# (which modules dominate at scale?)
# =============================================================================
max_n <- max(perproc_key$n_species)

stacked <- perproc_key |>
  filter(n_species == max_n) |>
  group_by(genome_size, phylogeny, quality, process) |>
  summarise(wall_time_min = mean(wall_time_min), .groups = "drop") |>
  mutate(
    run_label = paste0(genome_size, "\n", phylogeny, " / ", quality)
  )

# Colour palette for stages
n_stages   <- length(stage_order)
stage_cols <- setNames(
  colorRampPalette(c("#e6f2ff", "#084594"))(n_stages),
  stage_order
)

p_stack <- ggplot(stacked,
    aes(x = run_label, y = wall_time_min, fill = process)) +
  geom_col(width = 0.7) +
  facet_wrap(~genome_size, scales = "free_x") +
  scale_fill_manual(values = stage_cols, name = "Module") +
  labs(
    x = NULL,
    y = paste0("Wall time (min)  [n = ", max_n, " genomes]"),
    title = paste0("Time composition at n = ", max_n, " genomes")
  ) +
  theme_bench() +
  theme(
    axis.text.x  = element_text(size = 8, angle = 20, hjust = 1),
    legend.key.size = unit(0.4, "cm")
  )

ggsave(file.path(outdir, "fig4_time_composition.pdf"), p_stack,
       width = 10, height = 5, useDingbats = FALSE)
ggsave(file.path(outdir, "fig4_time_composition.png"), p_stack,
       width = 10, height = 5, dpi = 300)
message("Saved: fig4_time_composition")

message("\nAll figures written to: ", outdir)
