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
  # Genome acquisition
  "NCBIGENOMEDOWNLOAD",
  # Annotation preparation
  "AGAT_SPKEEPLONGESTISOFORM", "GFFREAD", "RENAME_FASTA",
  # Orthology
  "ORTHOFINDER_CAFE", "ORTHOFINDER_V2_CAFE",
  # Tree
  "RESCALE_TREE",
  # CAFE
  "CAFE_PREP", "CAFE_RUN_K", "CAFE_SELECT_K",
  "CAFE_RUN_BEST", "CAFE_MODEL_COMPARE", "CAFE_PLOT",
  # EggNOG / GO (only present when --run_eggnog used)
  "EGGNOG_DOWNLOAD", "EGGNOGMAPPER",
  "EGGNOG_TO_GO", "EGGNOG_TO_OG_GO", "OG_ANNOTATION_SUMMARY",
  "CAFE_GO_PREP", "CAFE_GO_RUN"
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

# Distinct, colourblind-friendly palette for pipeline stages
# (Okabe-Ito extended with a few extras)
stage_cols <- c(
  # Genome acquisition — blues
  "NCBIGENOMEDOWNLOAD"        = "#56B4E9",  # sky blue
  "GUNZIP"                    = "#0072B2",  # dark blue
  # Annotation preparation — greens
  "AGAT_SPKEEPLONGESTISOFORM" = "#009E73",  # teal
  "GFFREAD"                   = "#44AA99",  # green-teal
  "RENAME_FASTA"              = "#117733",  # dark green
  # Orthology — oranges
  "ORTHOFINDER_CAFE"          = "#E69F00",  # orange
  "ORTHOFINDER_V2_CAFE"       = "#D55E00",  # vermillion
  # Tree
  "RESCALE_TREE"              = "#999999",  # grey
  # CAFE — purples
  "CAFE_PREP"                 = "#CC79A7",  # pink
  "CAFE_RUN_K"                = "#882255",  # wine
  "CAFE_SELECT_K"             = "#AA4499",  # mauve
  "CAFE_RUN_BEST"             = "#332288",  # dark purple
  "CAFE_MODEL_COMPARE"        = "#6666AA",  # mid purple
  "CAFE_PLOT"                 = "#BBBBDD",  # pale purple
  # EggNOG / GO — yellows/browns
  "EGGNOG_DOWNLOAD"           = "#F0E442",  # yellow
  "EGGNOGMAPPER"              = "#DDCC77",  # sand
  "EGGNOG_TO_GO"              = "#AA8800",  # dark gold
  "EGGNOG_TO_OG_GO"           = "#886600",  # brown-gold
  "OG_ANNOTATION_SUMMARY"     = "#664400",  # dark brown
  "CAFE_GO_PREP"              = "#CC6677",  # rose
  "CAFE_GO_RUN"               = "#AA2233"   # dark red
)

# Label segments that are >= 5% of their bar total (so tiny slivers stay clean)
stacked_labelled <- stacked |>
  group_by(run_label, genome_size) |>
  mutate(
    bar_total  = sum(wall_time_min),
    pct        = wall_time_min / bar_total,
    seg_label  = ifelse(pct >= 0.05,
                        gsub("_", "\n", process),
                        NA_character_),
    # position label in the middle of its segment
    label_y    = cumsum(wall_time_min) - wall_time_min / 2
  ) |>
  ungroup()

p_stack <- ggplot(stacked_labelled,
    aes(x = run_label, y = wall_time_min, fill = process)) +
  geom_col(width = 0.7, colour = "white", linewidth = 0.3) +
  geom_text(aes(y = label_y, label = seg_label),
            size = 2.2, lineheight = 0.85, colour = "black",
            na.rm = TRUE) +
  facet_wrap(~genome_size, scales = "free_x") +
  scale_fill_manual(values = stage_cols, name = "Module",
                    drop = TRUE) +
  labs(
    x = NULL,
    y = paste0("Wall time (min)  [n = ", max_n, " genomes]"),
    title = paste0("Time composition at n = ", max_n, " genomes")
  ) +
  theme_bench() +
  theme(
    axis.text.x     = element_text(size = 8, angle = 20, hjust = 1),
    legend.key.size = unit(0.4, "cm"),
    legend.text     = element_text(size = 8)
  )

ggsave(file.path(outdir, "fig4_time_composition.pdf"), p_stack,
       width = 11, height = 6, useDingbats = FALSE)
ggsave(file.path(outdir, "fig4_time_composition.png"), p_stack,
       width = 11, height = 6, dpi = 300)
message("Saved: fig4_time_composition")

# =============================================================================
# Figure 5 — User guidance: will the pipeline work for my data?
#
# Shows the benchmark conditions in genome-size × assembly-quality space,
# coloured by wall time per genome and marked by CAFE convergence success.
# Reference zones for common taxa help users locate their own data.
# =============================================================================

# Approximate median genome sizes (Mb) for each benchmark category
genome_size_mb <- c(bacteria = 4, insect = 350, mammal = 3000)

# Approximate median scaffold N50 (kb) for each quality band per genome size
# (based on the N50 cutoffs used to build the input CSVs)
n50_kb_approx <- tribble(
  ~genome_size, ~quality,      ~n50_kb_approx,
  "bacteria",   "contiguous",  1000,
  "bacteria",   "fragmented",  40,
  "insect",     "contiguous",  5000,
  "insect",     "fragmented",  300,
  "mammal",     "contiguous",  120000,
  "mammal",     "fragmented",  15000
)

# Add convergence success flag — CAFE_RUN_BEST absent from per-process means
# CAFE did not converge for that run
converged_runs <- perproc |>
  group_by(run_id) |>
  summarise(cafe_converged = any(process == "CAFE_RUN_BEST"), .groups = "drop")

guidance <- metrics |>
  left_join(converged_runs, by = "run_id") |>
  left_join(n50_kb_approx, by = c("genome_size", "quality")) |>
  mutate(
    genome_mb       = genome_size_mb[as.character(genome_size)],
    min_per_genome  = total_wall_time_min / n_species,
    cafe_converged  = replace_na(cafe_converged, FALSE)
  )

# Reference zones: typical genome size and N50 ranges for common taxa
# (approximate — for orientation only)
ref_zones <- tribble(
  ~label,           ~xmin,  ~xmax,   ~ymin,    ~ymax,
  "Fungi",           8,      60,      10,       5000,
  "Plants\n(small)", 100,    500,     500,      50000,
  "Plants\n(large)", 1000,   20000,   50,       10000,
  "Fish",            400,    1500,    1000,     100000,
  "Birds",           1000,   2000,    10000,    500000
)

p_guidance <- ggplot() +
  # Reference zones (behind data)
  geom_rect(data = ref_zones,
            aes(xmin = xmin, xmax = xmax, ymin = ymin, ymax = ymax),
            fill = "grey90", colour = "grey70", linewidth = 0.4, alpha = 0.6,
            inherit.aes = FALSE) +
  geom_text(data = ref_zones,
            aes(x = sqrt(xmin * xmax), y = sqrt(ymin * ymax), label = label),
            size = 2.8, colour = "grey45", fontface = "italic", lineheight = 0.85,
            inherit.aes = FALSE) +
  # Benchmark data points
  geom_point(data = guidance,
             aes(x = genome_mb, y = n50_kb_approx,
                 colour = min_per_genome,
                 shape  = cafe_converged),
             size = 5, stroke = 1.2) +
  # Genome-size category labels
  geom_text(data = guidance |> distinct(genome_size, genome_mb),
            aes(x = genome_mb, y = 0.8, label = genome_size),
            size = 3, fontface = "bold", colour = "grey30",
            inherit.aes = FALSE) +
  scale_x_log10(
    name   = "Genome size (Mb)",
    labels = label_comma(),
    breaks = c(1, 4, 10, 100, 350, 1000, 3000, 10000)
  ) +
  scale_y_log10(
    name   = "Assembly scaffold N50 (kb)",
    labels = label_comma(),
    breaks = c(1, 10, 40, 100, 300, 1000, 5000, 15000, 50000, 120000)
  ) +
  scale_colour_gradient(
    low  = "#2166ac", high = "#d73027",
    name = "Wall time\nper genome (min)",
    labels = comma
  ) +
  scale_shape_manual(
    values = c("TRUE" = 16, "FALSE" = 4),
    name   = "CAFE converged",
    labels = c("TRUE" = "Yes", "FALSE" = "No")
  ) +
  facet_wrap(~phylogeny, labeller = label_both) +
  labs(
    title    = "Will the pipeline work for your data?",
    subtitle = paste0(
      "Benchmark points use representative median genome size and scaffold N50 for each category ",
      "(bacteria ~4 Mb, insect ~350 Mb, mammal ~3,000 Mb).\n",
      "These are approximate — locate your organism relative to the benchmark points ",
      "to estimate expected run time and likelihood of CAFE convergence."
    )
  ) +
  theme_bench() +
  theme(legend.box = "horizontal")

ggsave(file.path(outdir, "fig5_user_guidance.pdf"), p_guidance,
       width = 10, height = 6, useDingbats = FALSE)
ggsave(file.path(outdir, "fig5_user_guidance.png"), p_guidance,
       width = 10, height = 6, dpi = 300)
message("Saved: fig5_user_guidance")

message("\nAll figures written to: ", outdir)
