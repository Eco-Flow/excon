# EXCON Benchmarking Suite

Benchmarks CAFE performance across a 4-factor design for a methods paper.

## Factors

| Factor | Levels |
|--------|--------|
| Genome size | `bacteria` · `insect` · `mammal` |
| Phylogeny | `close` · `diverse` |
| Annotation quality | `contiguous` · `fragmented` |
| Dataset size | 10 · 50 · 100 species |

Each combination runs **CAFE only** — no GO enrichment, no genome quality stats.

---

## Quick start

The benchmark script runs on the **head/login node** exactly like any `nextflow run`
call — it submits jobs to the cluster via your scheduler. Just pass your usual
`--custom_config` and `-profile singularity`.

```bash
cd excon-1/

# HPC with SGE + Singularity (typical usage)
./benchmark/run_benchmark.sh \
    --profile singularity \
    --custom-config /path/to/your_hpc.config

# Subset for a quick sanity-check (bacteria n=10 only)
./benchmark/run_benchmark.sh \
    --profile singularity \
    --custom-config /path/to/your_hpc.config \
    --genome-sizes bacteria \
    --dataset-sizes 10

# Dry run — prints every Nextflow command without executing
./benchmark/run_benchmark.sh \
    --profile singularity \
    --custom-config /path/to/your_hpc.config \
    --dry-run

# Local Docker (laptop/workstation, no HPC)
./benchmark/run_benchmark.sh --profile docker
```

The `--custom-config` value is passed straight through as `--custom_config` to
every Nextflow invocation, so it works exactly the same as your normal runs.

Completed runs are recorded in `benchmark/results/run_log.tsv`. The script
skips any run already marked `COMPLETED`, so it is safe to re-run after a
crash — Nextflow `-resume` will also pick up from the last cached step within
each run.

---

## Options

| Flag | Default | Description |
|------|---------|-------------|
| `--genome-sizes` | `bacteria,insect,mammal` | Comma-separated subset |
| `--dataset-sizes` | `10,50,100` | N species to draw from the top of each CSV |
| `--phylogenies` | `close,diverse` | Subset |
| `--qualities` | `contiguous,fragmented` | Subset |
| `--profile` | `singularity` | Nextflow `-profile` (e.g. `singularity`, `docker`) |
| `--custom-config` | _(none)_ | Path to your HPC config file (sets executor, queues, cache dirs etc.) |
| `--max-memory` | `128.GB` | Memory cap passed to Nextflow |
| `--max-cpus` | `16` | CPU cap passed to Nextflow |
| `--no-resume` | off | Disable `-resume` (force fresh runs) |
| `--dry-run` | off | Print commands without running |

---

## Input CSVs

Each CSV in `inputs/` is named `{genome_size}_{phylogeny}_{quality}.csv`.
They contain 20 species each — **the runner takes the first N rows** for
each dataset-size level.

**Important:** Check that all accession numbers are current before running.
NCBI accessions occasionally get superseded. Verify with:

```bash
# Check one accession resolves
ncbi-genome-download --dry-run --assembly-accessions GCF_000006885.1 bacteria
```

### Species selection criteria

| Factor | Contiguous | Fragmented |
|--------|-----------|------------|
| Bacteria | Complete chromosome assemblies; typically 1–3 sequences, N50 ≥ 500 kb | WGS assemblies; >50 scaffolds, N50 < 100 kb |
| Insect | Chromosome-level (scaffold N50 > 1 Mb) | Contig/scaffold assemblies (N50 < 100 kb) |
| Mammal | Chromosome-level RefSeq reference genomes | Contig or scaffold assemblies |

For 50 and 100 species you will need to add more rows to each CSV —
just append `species_name,GCF_XXXXXXXXX.X` lines, preserving the ordering
(contiguous = better assemblies at the top, fragmented = worse at the top).

---

## Outputs

After each run:
```
benchmark/results/
├── run_log.tsv                        # One row per run, completion status
├── {run_id}.log                       # Raw Nextflow stdout/stderr
├── {run_id}_input.csv                 # Subset CSV used for that run
├── {run_id}/                          # Nextflow outdir
│   └── pipeline_info/
│       ├── execution_trace.tsv        # Per-task metrics (parsed below)
│       ├── execution_timeline.html    # Gantt chart
│       └── execution_report.html      # Resource summary
└── {run_id}_work/                     # Nextflow work dir (for -resume)
```

After all runs complete, the runner calls `collect_metrics.py` automatically.

### Metrics files

| File | Contents |
|------|---------|
| `benchmark_metrics.tsv` | One row per run: wall time, CPU-hours, peak RAM, efficiency, throughput |
| `benchmark_per_process.tsv` | One row per (run × process): per-stage breakdown |

### Figures produced by plot_benchmark.R

| File | Description |
|------|-------------|
| `fig1_overall_scaling.pdf` | Wall time + peak RAM vs n_species; lines per genome size, faceted by phylogeny |
| `fig2_per_module_scaling.pdf` | One panel per pipeline stage; wall time vs n_species |
| `fig3_efficiency.pdf` | CPU efficiency % + genomes/hour vs n_species |
| `fig4_time_composition.pdf` | Stacked bar of module time at the largest N |

R packages required: `ggplot2`, `dplyr`, `tidyr`, `scales`, `patchwork`

### Metric definitions

| Metric | Definition |
|--------|-----------|
| `total_wall_time_min` | Clock time from first task submit to last task completion |
| `total_cpu_hours` | Σ (realtime × %CPU/100) across all tasks |
| `peak_ram_gb` | Maximum `peak_rss` across all tasks in the run |
| `mean_cpu_pct` | Realtime-weighted mean CPU utilisation (%) |
| `cpu_efficiency` | `cpu_hours / (max_cpus × wall_hours)` |
| `genomes_per_hour` | `n_species / wall_hours` |

---

## Workflow: run now, collect/plot later

The runner, metrics collector, and plotter are fully independent — you can
run benchmark batches whenever you like and collect/plot only when you have
enough results to compare.

```bash
# 1. Run some benchmarks (can be called multiple times, adds to run_log.tsv)
./benchmark/run_benchmark.sh --profile singularity --custom-config my_hpc.config \
    --genome-sizes bacteria --dataset-sizes 10,50

./benchmark/run_benchmark.sh --profile singularity --custom-config my_hpc.config \
    --genome-sizes insect --dataset-sizes 10,50

# 2. Collect metrics from all COMPLETED runs so far
python3 benchmark/collect_metrics.py \
    --results-dir benchmark/results \
    --log         benchmark/results/run_log.tsv \
    --output      benchmark/results/benchmark_metrics.tsv \
    --max-cpus    16

# 3. Plot (requires R packages: ggplot2, dplyr, tidyr, scales, patchwork)
Rscript benchmark/plot_benchmark.R \
    --metrics benchmark/results/benchmark_metrics.tsv \
    --perproc benchmark/results/benchmark_per_process.tsv \
    --outdir  benchmark/results/figures
```

---

## Computational scale notes

| Genome size | ~Genome size | 10 species | 50 species | 100 species |
|-------------|-------------|-----------|-----------|------------|
| Bacteria | 0.5–10 Mb | ~15 min | ~1 h | ~3–5 h |
| Insect | 150–2000 Mb | ~2 h | ~12–24 h | ~days |
| Mammal | 2000–3500 Mb | ~8 h | ~days | very long |

OrthoFinder scales super-linearly with species count and genome size.
For mammals at 50–100 species, consider using `--input_tree` and
`--input_orthogroups` to reuse a pre-computed OrthoFinder run.
