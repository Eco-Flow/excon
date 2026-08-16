# CAFE Gene Family Evolution Results — Data Guide

This archive contains the outputs from the EXCON pipeline, which runs CAFE5 to detect gene family
expansion and contraction across a set of species, and follows up with GO
enrichment analysis on the significantly evolving families.

Pipeline version: EXCON v2.4.0  
CAFE5 version: 4.2.1  
OrthoFinder version: see `pipeline_info/software_versions.yml`

---

## Directory overview

```
cafe/                   CAFE5 model fitting and selection
cafe_plot/              Tree visualisations of expansion/contraction
cafe_go/                Per-species GO enrichment on CAFE results
cafe_go_summary/        Cross-species GO summary heatmaps and dot plots
cafe_go_large/          GO enrichment for high-differential families (if run)
```

---

## `cafe/`

Contains the raw CAFE5 outputs organised by model.

### `cafe/base/`

Results from the preparatory CAFE5 run used to estimate the global gene-family
turnover rate (lambda) and fit an error model for annotation noise.

| File | Description |
|------|-------------|
| `hog_gene_counts.tsv` | Input gene-family count table passed to CAFE5. One row per orthogroup/HOG, one column per species, values are gene copy numbers. |
| `SpeciesTree_rooted_ultra.txt` | Ultrametric rooted species tree (Newick) used by CAFE5. Branch lengths are in units scaled for CAFE5 (see `--tree_scale_factor` parameter). |
| `lambda.txt` | Global lambda estimate from the base CAFE5 run. Lambda is the per-gene per-unit-time rate of gene gain and loss assumed under the birth-death model. |
| `hog_filtering_report.tsv` | (optional) Summary of gene families filtered before CAFE5, if `--max_differential` filtering was applied. Columns: family ID, max pairwise species difference, kept/removed status. |
| `hog_gene_counts_large.tsv` | (optional) Families excluded from the main run because their copy-number variance was too high for one shared model to fit. Each is analysed on its own (see `cafe/large_families/`). |
| `Out_cafe/Base_count.tab` | Gene counts at each node inferred by CAFE5 under the base (uniform-lambda) model during error-model fitting. |
| `Out_cafe_errormodel/Base_error_model.txt` | Fitted error model parameters. CAFE5 uses this to account for gene-count errors caused by incomplete genome assemblies or annotation artefacts. |
| `cafe_base.log` | CAFE5 stdout/stderr log from the baseline run. |
| `cafe_errormodel.log` | CAFE5 stdout/stderr log from the error-model fitting run. |
| `pruned_tree/` | Intermediate tree files produced during species-tree processing. |

### `cafe/gamma/` and `cafe/gamma_per_family/`

CAFE5 model runs with rate variation. CAFE5 fits two model variants and the
pipeline selects the better one by AIC (see `cafe/model_comparison/`).

| File | Description |
|------|-------------|
| `Gamma_results.txt` | Overall model fit summary: log-likelihood, lambda, alpha (gamma shape). |
| `Gamma_family_results.txt` | Per-family p-values for the hypothesis that the family evolved under the fitted rate (significant families have p ≤ 0.05 by default). |
| `Gamma_branch_probabilities.tab` | Matrix of branch-level p-values. Rows = gene families, columns = species/nodes. Values are the probability that the observed change on that branch occurred by chance under the null model. Used to identify which branch a family significantly changed on. |
| `Gamma_change.tab` | Matrix of inferred gene count changes. Rows = gene families, columns = species/nodes. Positive = genes gained on that branch, negative = genes lost. These changes are ancestral state reconstructions, not direct observations. |
| `Gamma_count.tab` | Absolute inferred gene counts at each node (ancestral + extant). |
| `Gamma_asr.tre` | Annotated phylogenetic tree in Newick format with inferred gene counts embedded. Used by cafeplotter and `cafe_plot_altviz.R`. |

The `gamma_per_family` subdirectory has the same files with per-family rate
variation (`-p` flag in CAFE5), where each family can have its own lambda.

### `cafe/model_comparison/`

| File | Description |
|------|-------------|
| `cafe_model_comparison.tsv` | AIC scores for each model tested (uniform lambda, gamma, gamma per-family). Lower AIC = better fit. |
| `best_model.txt` | Name of the selected model (`gamma` or `gamma_per_family`). |

### `cafe/best/`

A copy of the CAFE5 output directory from the best-fitting model. This is the
canonical result used for all downstream analysis (visualisation, GO
enrichment). Contents are identical to either `cafe/gamma/` or
`cafe/gamma_per_family/` depending on model selection.

### `cafe/large_families/`

(optional) CAFE5 results for gene families that were excluded from the main run
because of extreme copy-number variation between species — one shared lambda
cannot fit both a modest and an extreme size differential at once, so these
are run independently, each fitting its own lambda:

| File/dir | Description |
|------|-------------|
| `per_family/Out_cafe_large_<HOG>/` | Raw single-family CAFE5 output for each high-differential HOG, run on its own. |
| `Out_cafe_large/` | All converged per-family runs stitched back into one CAFE5-shaped directory (same file structure as `cafe/gamma/`), so it can be read like an ordinary CAFE5 result by `cafe_plot_large/`, `cafe_go_large/`, etc. |
| `large_family_lambda_summary.tsv` | Each family's own independently-fitted lambda and -lnL, so it is visible which rate was used where instead of one shared value. |

Families that still fail to converge even on their own (rare) are dropped from
`Out_cafe_large/` rather than causing the merge to fail.

---

## `cafe_plot/`

Tree visualisations of the expansion/contraction results. All figures are
provided in both PDF and SVG format. SVG files are suitable for editing in
Inkscape or Illustrator.

### Cafeplotter summary tree (`cafe_plotter/`)

Produced by [CAFEPlotter](https://github.com/moshi4/CafePlotter).

| File | Description |
|------|-------------|
| `cafe_plotter/summary_all_gene_family.pdf` / `.svg` | Species tree annotated with the total number of gene families expanding (+) and contracting (−) on each branch. **No p-value filter is applied** — all families with any directional change are counted. This is the standard CAFE summary visualisation. |
| `cafe_plotter/gene_family/` | One plot per significantly evolving gene family (family-level p ≤ 0.05). Each shows the inferred gene count at every node, with asterisks for branch-level significance and red/blue colouring for expansion/contraction. |

### Alternative summary trees

Two additional tree figures are produced by the pipeline using data from
`CAFE_summary.txt` (see `cafe_go/`). Unlike the cafeplotter output, these are
restricted to **statistically significant branches** (branch p ≤ 0.05).

| File | Description |
|------|-------------|
| `cafe_sig_hog_tree.pdf` / `.svg` | Species tree annotated with the number of **significantly evolving gene families** per branch (+expanded / −contracted, p ≤ 0.05). Use this figure when reporting how many families changed on each lineage at the chosen significance threshold. |
| `cafe_sig_gene_tree.pdf` / `.svg` | Species tree annotated with the **net number of genes gained or lost** in significant families per branch. A family that shrank from 10 to 7 genes contributes −3. Use this to compare the magnitude of gene gain/loss across lineages, not just family counts. |

### Node number guide

| File | Description |
|------|-------------|
| `cafe_node_label_guide.pdf` / `.svg` | Species tree with internal node numbers labelled. CAFE5 assigns numbers to ancestral nodes (e.g., node 14, 15, 16), and these numbers appear throughout the results tables. Use this figure as a key to identify which ancestral lineage each node number refers to. |

---

## `cafe_go/`

Per-species GO enrichment analysis on gene families identified by CAFE5. The
subdirectory name encodes the parameters used:

```
cafe_go/{algorithm}_cutoff{p}_type{ontology}/
```

For example: `cafe_go/classic_fisher_cutoff0.05_typeBP_MF_CC/`

Parameters:
- `algorithm`: topGO algorithm used (`classic_fisher` or `weight01`; see topGO documentation for the distinction)
- `cutoff`: p-value threshold applied to CAFE5 branch probabilities to define significant families
- `type`: GO ontology categories tested (BP = Biological Process, MF = Molecular Function, CC = Cellular Component)

### `CAFE_summary.txt`

A tab-separated summary table — the key quantitative result of the CAFE5
analysis. One row per species (leaf node) and per internal ancestral node.
Comment lines at the top of the file explain every column.

| Column | Description |
|--------|-------------|
| `Species/Node` | Leaf species name or internal node number (node numbers match `cafe_node_label_guide.pdf`) |
| `Total_HOGs_significant` | Total gene families with a significant branch p-value |
| `Expansion_HOGs_significant` | Significant families with a net positive copy-number change |
| `Contraction_HOGs_significant` | Significant families with a net negative copy-number change |
| `Expansion_genes_significant` | Sum of copy-number increases across all significant expanding families (total genes gained) |
| `Contraction_genes_significant` | Sum of copy-number decreases across all significant contracting families (reported as a negative number; total genes lost) |
| `Expansion_HOGs_total` | All families with any expansion regardless of significance — corresponds to the (+) counts shown in the cafeplotter summary tree |
| `Contraction_HOGs_total` | All families with any contraction regardless of significance — corresponds to the (−) counts in the cafeplotter summary tree |

Note on terminology: "HOG" is used throughout as a generic label for gene
family. When OrthoFinder v2 is used (`--orthofinder_v2`), these are
Hierarchical Orthogroups (HOGs) from `N0.tsv`. When OrthoFinder v3 is used
(default), these are flat Orthogroups from `Orthogroups.tsv`. Both are handled
identically by the pipeline.

### `go_inputs/`

Raw input files used for the topGO enrichment runs. Provided for
reproducibility — these files can be used to re-run or modify the GO analysis
independently of the pipeline.

| File pattern | Description |
|---|---|
| `{species}.pos.txt` | Orthogroup/HOG IDs for significantly **expanding** families in this species (target list for GO enrichment). |
| `{species}.neg.txt` | Orthogroup/HOG IDs for significantly **contracting** families in this species. |
| `{species}.BK.txt.uniq` | Background set of all orthogroup/HOG IDs detected in this species (universe for Fisher's exact test). Unique entries only. |

### GO enrichment results (per species)

One set of files per species per direction (positive = expanding, negative =
contracting). Files are named `{species}.{pos|neg}_TopGo_results_ALL.tab` and
associated plots.

| File pattern | Description |
|---|---|
| `{species}.{pos\|neg}_TopGo_results_ALL.tab` | Full topGO results table. Columns include GO ID, term name, number of annotated genes in the universe, number of significant genes, p-value, Bonferroni-corrected p-value, and fold enrichment. All GO terms tested are included regardless of significance. |
| `TopGO_barplot_{species}.{pos\|neg}.pdf` / `.svg` | Horizontal bar chart of the top significant GO terms, grouped by ontology (BP/MF/CC), coloured by category, x-axis shows −log₁₀(p-value). |
| `TopGO_dotplot_{species}.{pos\|neg}.pdf` / `.svg` | Dot plot of the same terms: x-axis = fold enrichment, dot size = number of significant genes annotated to the term, dot colour gradient = −log₁₀(p-value). |

---

## `cafe_go_summary/`

Cross-species GO enrichment summary visualisations. These figures combine the
per-species results from `cafe_go/` into a single comparative view.

### `cafe_go_summary/cafe_go/`

Summary for the main CAFE5 run.

| File | Description |
|------|-------------|
| `go_enrichment_heatmap.pdf` / `.svg` / `.png` | Heatmap of GO term enrichment across all species, with terms on rows and species on columns. Colour intensity encodes −log₁₀(p-value). Expanding and contracting families are shown in separate panels. |
| `go_enrichment_dotplot.pdf` / `.svg` / `.png` | Dot plot version of the above: dot size = number of significant genes, colour = significance. |
| `go_enrichment_heatmap_aligned.pdf` / `.svg` / `.png` | As above but with species and terms ordered to align shared patterns. |
| `go_enrichment_heatmap_aligned_shared.pdf` / `.svg` / `.png` | Shows only GO terms significant in two or more species. Useful for identifying shared evolutionary signals. |
| `go_enrichment_heatmap_expanding_only.pdf` / `.svg` / `.png` | Summary heatmap restricted to expanding families. |
| `go_enrichment_heatmap_contracting_only.pdf` / `.svg` / `.png` | Summary heatmap restricted to contracting families. |
| `go_enrichment_dotplot_aligned.pdf` / `.svg` / `.png` | Aligned dot plot across all species. |

| File | Description |
|------|-------------|
| `Go_summary_pos.pdf` / `.svg` | Heatmap across species for GO terms enriched in **expanding** families (rows = GO terms, columns = species). |
| `Go_summary_neg.pdf` / `.svg` | Equivalent for **contracting** families. |
| `Go_summary_pos_noNode.pdf` / `.svg` | As above but with ancestral nodes excluded — shows leaf species only. |
| `Go_summary_neg_noNode.pdf` / `.svg` | Contracting, leaf species only. |

### `cafe_go_summary/cafe_go/raw_plotting_files/`

Data tables and R scripts used to produce the summary figures. Provided for
reproducibility and to allow custom re-plotting.

| File | Description |
|------|-------------|
| `Go_summary_pos.tsv` | Combined GO enrichment table for expanding families across all species. Input to the summary heatmap plots. |
| `Go_summary_neg.tsv` | Equivalent for contracting families. |
| `Go_summary_posneg_merged.tsv` | Merged table of both directions. |
| `plotting_go.R` | R script that reads the above TSV files and generates the heatmap / dot-plot summary figures. |
| `plotting_go_summary.R` | R script that generates the multi-panel enrichment summary figures. |

### `cafe_go_summary/cafe_go_large/`

(optional) Same structure as `cafe_go_summary/cafe_go/` but for the
high-differential family subset (see `cafe/large_families/`).

---

## `cafe_go_large/`

(optional) Per-species GO enrichment on the high-differential gene families
that were excluded from the main CAFE5 run. Same structure as `cafe_go/`.

---

## Reproducing the analysis

To regenerate all outputs from the raw genome annotations, run:

```bash
nextflow run main.nf \
  --input samplesheet.csv \
  --outdir results \
  --predownloaded_gofiles /path/to/go_files \
  -profile singularity
```

See the pipeline README and `docs/outputs.md` for full parameter descriptions
and output documentation.

The exact software versions used are recorded in
`pipeline_info/software_versions.yml`.
