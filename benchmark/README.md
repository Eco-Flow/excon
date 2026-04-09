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

### Running on HPC

The benchmark script calls `nextflow run` sequentially for each factor
combination — the full suite could run for days, so it should not be left in
a plain terminal session. Two options:

**Option A — submit as a job (recommended)**

The benchmark script itself uses almost no resources (it just submits child
jobs via your scheduler). Give it a long walltime and minimal memory:

```bash
# submit_benchmark.sh
#!/usr/bin/env bash
#$ -N excon_benchmark
#$ -l h_rt=240:00:00   # long walltime — full suite can take days
#$ -l h_vmem=4G        # low memory — just coordinates job submission
#$ -pe smp 1
#$ -cwd
#$ -j y
#$ -o benchmark_runner.log

./benchmark/run_benchmark.sh \
    --profile singularity \
    --custom-config /path/to/your_hpc.config
```

```bash
qsub submit_benchmark.sh
```

Add `--parallel` to launch all runs simultaneously rather than one after the
other — on HPC this is strongly recommended, as each run submits its own jobs
to the scheduler independently:

```bash
#$ ...
./benchmark/run_benchmark.sh \
    --profile singularity \
    --custom-config /path/to/your_hpc.config \
    --parallel
```

Nextflow running inside a compute job can still call `qsub` to submit its own
child jobs — this is standard practice and schedulers allow it. The process
hierarchy looks like:

```
qsub submit_benchmark.sh        ← low-resource, long-walltime job
  └─ nextflow run main.nf        ← Nextflow driver (lightweight)
       ├─ qsub NCBIGENOMEDOWNLOAD
       ├─ qsub ORTHOFINDER
       ├─ qsub CAFE_PREP
       └─ ...                    ← actual compute on worker nodes
```

**Option B — screen/tmux on the head node**

If your HPC allows long-running head node processes (Nextflow itself is
lightweight, so many sysadmins permit this):

```bash
screen -S benchmark
./benchmark/run_benchmark.sh \
    --profile singularity \
    --custom-config /path/to/your_hpc.config
# Ctrl+A D to detach; screen -r benchmark to reattach
```

Either way, because completed runs are logged in `run_log.tsv` and Nextflow
uses `-resume`, any interruption can be recovered by simply re-running the
same command.

### Subset runs and dry run

```bash
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

# Local Docker on Mac (requires --orthofinder_v2 and mac profile)
./benchmark/run_benchmark.sh \
    --profile docker,mac \
    --nf-args "--orthofinder_v2"
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
| `--parallel` | off | Launch all runs simultaneously (recommended on HPC) |
| `--no-resume` | off | Disable `-resume` (force fresh runs) |
| `--dry-run` | off | Print commands without running |
| `--nf-args` | _(none)_ | Extra pipeline flags passed through verbatim (quoted), e.g. `"--orthofinder_v2"` |

---

## Input CSVs

Each CSV in `inputs/` is named `{genome_size}_{phylogeny}_{quality}.csv`.
They contain 20 species each — **the runner takes the first N rows** for
each dataset-size level, so ordering matters:

- **contiguous** files: best-quality assemblies at the top
- **fragmented** files: most-fragmented assemblies at the top

For n=50 and n=100 you will need to add more rows. See [species selection criteria](#species-selection-criteria) below.

### Validating accessions before running

The provided CSVs are a starting point but **some accessions may be wrong or
withdrawn** — NCBI accessions can be superseded, and version numbers (`.1`,
`.2`) change over time. Run the validator before your first benchmark run:

```bash
# Requires 'datasets' (NCBI CLI) or 'ncbi-genome-download' in your environment
./benchmark/validate_inputs.sh

# Or check one file at a time
./benchmark/validate_inputs.sh benchmark/inputs/bacteria_close_contiguous.csv
```

This prints `OK` / `FAIL` for every accession and lists all broken ones at the
end. Fix any `FAIL` entries before running — a bad accession causes the whole
Nextflow run to abort with:

```
ERROR: No downloads matched your filter. Please check your options.
```

To find the correct accession for a species, use the NCBI datasets CLI:

```bash
datasets summary genome taxon "Streptococcus pyogenes" \
    --assembly-level chromosome \
    --as-json-lines | \
    dataformat tsv genome --fields accession,organism-name
```

### Recovering from a bad accession mid-run

If a run fails due to a bad accession:

1. Fix the accession in the relevant CSV
2. Remove the `FAILED` entry from `benchmark/results/run_log.tsv` so the runner retries it
3. Rerun — Nextflow `-resume` will skip already-completed tasks

```bash
# Remove failed entries from the log
grep -v "FAILED" benchmark/results/run_log.tsv > /tmp/log_clean.tsv \
    && mv /tmp/log_clean.tsv benchmark/results/run_log.tsv
```

### Nextflow session lock errors

If you see `Unable to acquire lock on session` when starting a run:

```bash
# Check if a previous Nextflow process is still running
lsof /Users/cwyatt/Downloads/excon-1/.nextflow/cache/*/db/LOCK

# If a Java PID is listed, kill it (replace 12345 with actual PID)
kill 12345
```

Each benchmark run uses its own isolated cache directory
(`benchmark/results/{run_id}_cache/`) so runs can never lock each other out.
The lock error only occurs if a separate Nextflow process from outside the
benchmark script is still active.

### Species selection criteria and N50 cutoffs

N50 is used to define the quality band for each CSV. The `--min-n50` and
`--max-n50` flags in `ncbi_table_to_csv.py` filter by the Scaffold N50 column
from the NCBI Assembly table.

| Factor | Contiguous (`--min-n50`) | Fragmented (`--min-n50` / `--max-n50`) |
|--------|--------------------------|----------------------------------------|
| Bacteria | `--min-n50 500000` (≥500 kb) | `--min-n50 10000 --max-n50 100000` |
| Insect | `--min-n50 2000000` (≥1 Mb) | `--min-n50 1000 --max-n50 1500000` |
| Mammal | `--min-n50 90000000` (≥90 Mb) | `--min-n50 3000000 --max-n50 50000000` |

The floor on the fragmented band (`--min-n50`) is important — without it you
include completely broken assemblies (thousands of tiny contigs) that cause
OrthoFinder or CAFE to fail.

### How the current CSVs were generated

The commands below document exactly how each input file was built, so the
species selection is reproducible. In each case the source TSV was downloaded
from the NCBI Assembly web interface (**Send to → File → ID Table**).


Mammals, close:Primates, diverse:Mammals
Insects, close:Drosophilidae, diverse:Diptera
Bacteria, close:Streptococcus, diverse:Alphaproteobacteria

e.g.:

**Mammal — diverse — contiguous** (`mammal_diverse_contiguous.csv`)
```bash
# Source: all Mammalia genomes from NCBI Assembly
python3 benchmark/ncbi_table_to_csv.py "ncbi_dataset (5).tsv" \
    --min-n50 90000000 \
    --output benchmark/inputs/mammal_diverse_contiguous.csv
```

**Mammal — diverse — fragmented** (`mammal_diverse_fragmented.csv`)
```bash
# Source: all Mammalia genomes from NCBI Assembly
python3 benchmark/ncbi_table_to_csv.py "ncbi_dataset (5).tsv" \
    --min-n50 3000000 --max-n50 50000000 \
    --output benchmark/inputs/mammal_diverse_fragmented.csv
```

**Primates — close — contiguous** (`mammal_close_contiguous.csv`)
```bash
# Source: all Mammalia genomes from NCBI Assembly
python3 benchmark/ncbi_table_to_csv.py "ncbi_dataset (4).tsv" \
    --min-n50 90000000 \
    --output benchmark/inputs/mammal_close_contiguous.csv
```

**Primates — close — fragmented** (`mammal_close_fragmented.csv`)
```bash
# Source: all Mammalia genomes from NCBI Assembly
python3 benchmark/ncbi_table_to_csv.py "ncbi_dataset (4).tsv" \
    --min-n50 3000000 --max-n50 50000000 \
    --output benchmark/inputs/mammal_close_fragmented.csv
```


### Building or extending a species list

1. Go to [https://www.ncbi.nlm.nih.gov/assembly](https://www.ncbi.nlm.nih.gov/assembly)
2. Search for your taxonomic group, e.g.:
   - `Streptococcus[Organism]` — all Streptococcus assemblies
   - `Mammalia[Organism]` — all mammals
   - `Proteobacteria[Organism]` — for diverse bacteria (combine multiple phyla)
3. Click **Send to** → **File** → Format: **ID Table (text)** → **Create File**
4. Convert with N50 filtering:

```bash
# Contiguous — one representative per species, high N50
python3 benchmark/ncbi_table_to_csv.py ncbi_table.tsv \
    --min-n50 500000 \
    --output benchmark/inputs/bacteria_close_contiguous.csv

# Fragmented — quality band with floor and ceiling
python3 benchmark/ncbi_table_to_csv.py ncbi_table.tsv \
    --min-n50 10000 --max-n50 100000 \
    --output benchmark/inputs/bacteria_close_fragmented.csv

# Append more species to an existing file
python3 benchmark/ncbi_table_to_csv.py more_taxa.tsv \
    --min-n50 500000 \
    --append --output benchmark/inputs/bacteria_diverse_contiguous.csv
```

The script keeps only `GCF_` (RefSeq) entries, picks the highest accession
number per species binomial (proxy for most recent), removes duplicates, and
converts organism names to underscore format.

Then validate before running:

```bash
./benchmark/validate_inputs.sh benchmark/inputs/bacteria_close_contiguous.csv
```

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
