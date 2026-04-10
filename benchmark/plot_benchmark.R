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

metrics_file  <- get_arg("--metrics",  "benchmark/results/benchmark_metrics.tsv")
perproc_file  <- get_arg("--perproc",  "benchmark/results/benchmark_per_process.tsv")
outdir        <- get_arg("--outdir",   "benchmark/results/figures")
metadata_file <- get_arg("--metadata", "benchmark/inputs/metadata.tsv")

if (!file.exists(metrics_file)) stop("Cannot find --metrics file: ", metrics_file)
if (!file.exists(perproc_file)) stop("Cannot find --perproc file: ", perproc_file)
dir.create(outdir, recursive = TRUE, showWarnings = FALSE)

# ---- Load metadata ----------------------------------------------------------
# metadata.tsv maps each (genome_size, phylogeny, quality) combination to a
# specific clade name and representative genome size / N50 for Fig 5.
# Edit benchmark/inputs/metadata.tsv to update these values.
if (file.exists(metadata_file)) {
  metadata <- read.delim(metadata_file, stringsAsFactors = FALSE)
} else {
  warning("metadata.tsv not found — clade labels and genome/N50 values will ",
          "fall back to generic category names. Expected: ", metadata_file)
  metadata <- data.frame(
    genome_size = character(), phylogeny = character(), quality = character(),
    clade = character(), genome_mb = numeric(), n50_kb = numeric()
  )
}

# ---- Load data --------------------------------------------------------------
## Genome size display order — add new categories here to control plot order
genome_size_levels <- c("bacteria", "insect", "mammal", "bird", "reptile",
                        "fish", "plant", "fungi")

metrics <- read.delim(metrics_file, stringsAsFactors = FALSE) |>
  mutate(
    genome_size = factor(genome_size,
                         levels = intersect(genome_size_levels,
                                            unique(genome_size))),
    quality     = factor(quality,     levels = c("contiguous", "fragmented")),
    phylogeny   = factor(phylogeny,   levels = c("close", "diverse")),
    n_species   = as.integer(n_species)
  ) |>
  left_join(metadata, by = c("genome_size", "quality", "phylogeny"))

perproc <- read.delim(perproc_file, stringsAsFactors = FALSE) |>
  mutate(
    genome_size = factor(genome_size,
                         levels = intersect(genome_size_levels,
                                            unique(genome_size))),
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

## Colours per genome-size category — add a new entry here for each new group.
## Any category not listed gets an auto-assigned grey as a fallback.
genome_colours_defined <- c(
  bacteria = "#4DAF4A",
  insect   = "#377EB8",
  mammal   = "#E41A1C",
  bird     = "#FF7F00",
  reptile  = "#984EA3",
  fish     = "#A65628",
  plant    = "#00CED1",
  fungi    = "#FFD700"
)
# Keep only categories present in the data; add grey for any unknown category
genome_size_present <- unique(as.character(metrics$genome_size))
auto_grey <- setdiff(genome_size_present, names(genome_colours_defined))
genome_colours <- c(
  genome_colours_defined[intersect(names(genome_colours_defined),
                                   genome_size_present)],
  setNames(rep("grey50", length(auto_grey)), auto_grey)
)
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
# Use each condition's own maximum available N, not a single global max.
# This ensures all conditions appear even when dataset sizes are uneven.
stacked <- perproc_key |>
  group_by(genome_size, phylogeny, quality) |>
  filter(n_species == max(n_species)) |>
  ungroup() |>
  group_by(genome_size, phylogeny, quality, process) |>
  summarise(
    wall_time_min = mean(wall_time_min),
    n_used        = first(n_species),
    .groups = "drop"
  ) |>
  mutate(
    run_label = paste0(genome_size, "\n", phylogeny, " / ", quality,
                       "\n(n=", n_used, ")")
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
    x        = NULL,
    y        = "Wall time (min)",
    title    = "Time composition at maximum available dataset size",
    subtitle = "n shown in each bar label — varies by category depending on available species"
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

# Convergence flag — CAFE_RUN_BEST absent means CAFE did not converge
converged_runs <- perproc |>
  group_by(run_id) |>
  summarise(cafe_converged = any(process == "CAFE_RUN_BEST"), .groups = "drop")

# metrics already has genome_mb, n50_kb, clade joined from metadata above
guidance <- metrics |>
  left_join(converged_runs, by = "run_id") |>
  mutate(
    # Fall back to generic label if metadata not supplied
    clade          = ifelse(is.na(clade) | clade == "",
                            paste(genome_size, phylogeny, sep = "\n"), clade),
    genome_mb      = ifelse(is.na(genome_mb), NA_real_, genome_mb),
    n50_kb         = ifelse(is.na(n50_kb),
                            c(contiguous = 1000, fragmented = 100)[
                              as.character(quality)], n50_kb),
    min_per_genome = total_wall_time_min / n_species,
    cafe_converged = replace_na(cafe_converged, FALSE)
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
             aes(x = genome_mb, y = n50_kb,
                 colour = min_per_genome,
                 shape  = cafe_converged),
             size = 5, stroke = 1.2) +
  # Clade labels next to each point
  geom_text(data = guidance |> distinct(genome_mb, n50_kb, clade),
            aes(x = genome_mb, y = n50_kb, label = clade),
            size = 2.5, hjust = -0.15, vjust = 0.5,
            lineheight = 0.85, colour = "grey20",
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
      "Genome size and scaffold N50 are representative values defined in ",
      "benchmark/inputs/metadata.tsv — not measured from assemblies.\n",
      "Locate your organism relative to the labelled benchmark clades to ",
      "estimate run time and likelihood of CAFE convergence."
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
