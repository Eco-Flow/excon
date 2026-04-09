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

### Step 1 — create run directories

`run_benchmark.sh` prepares a self-contained directory for each factor
combination. It does not run Nextflow itself.

```bash
# HPC — creates run dirs with your profile and config baked in
./benchmark/run_benchmark.sh \
    --profile singularity \
    --custom-config /path/to/your_hpc.config

# Subset — bacteria n=10 only
./benchmark/run_benchmark.sh \
    --profile singularity \
    --custom-config /path/to/your_hpc.config \
    --genome-sizes bacteria \
    --dataset-sizes 10

# Local Mac with Docker
./benchmark/run_benchmark.sh \
    --profile docker,mac \
    --nf-args "--orthofinder_v2"
```

Each combination produces `benchmark/results/{run_id}/`:
```
benchmark/results/bacteria_close_contiguous_n10/
    input.csv   ← 10-species samplesheet
    run.sh      ← nextflow command, ready to execute
```

### Step 2 — run

Nextflow must run on the **login node** (so it can submit child jobs to the
scheduler). Run `run.sh` however you like — each run is fully isolated so you
can run multiple simultaneously.

```bash
# One run interactively
cd benchmark/results/bacteria_close_contiguous_n10
./run.sh

# Resume after interruption
./run.sh -resume

# Two runs in parallel — nohup keeps them alive after you disconnect
nohup benchmark/results/bacteria_close_contiguous_n10/run.sh \
    > benchmark/results/bacteria_close_contiguous_n10/nextflow.log 2>&1 &

nohup benchmark/results/bacteria_close_fragmented_n10/run.sh \
    > benchmark/results/bacteria_close_fragmented_n10/nextflow.log 2>&1 &

# Follow progress
tail -f benchmark/results/bacteria_close_contiguous_n10/nextflow.log
qstat   # see child jobs on the scheduler
```

The process hierarchy for each run:
```
login node: nohup run.sh             ← lightweight, stays on login node
  └─ nextflow run main.nf            ← Nextflow driver (lightweight)
       ├─ qsub NCBIGENOMEDOWNLOAD    ← actual compute on worker nodes
       ├─ qsub ORTHOFINDER
       └─ ...
```

### Step 3 — collect metrics

Once any runs have completed:

```bash
python3 benchmark/collect_metrics.py \
    --results-dir benchmark/results \
    --output      benchmark/results/benchmark_metrics.tsv
```

The collector scans for completed runs automatically (any directory containing
`output/pipeline_info/execution_trace.tsv`). Run it at any point — it picks up
whatever has finished so far.

---

## Options

| Flag | Default | Description |
|------|---------|-------------|
| `--genome-sizes` | `bacteria,insect,mammal` | Comma-separated subset |
| `--dataset-sizes` | `10,50,100` | N species to draw from the top of each CSV |
| `--phylogenies` | `close,diverse` | Subset |
| `--qualities` | `contiguous,fragmented` | Subset |
| `--profile` | `singularity` | Nextflow `-profile` (e.g. `singularity`, `docker`) |
| `--custom-config` | _(none)_ | Path to your HPC config file |
| `--max-memory` | `128.GB` | Memory cap passed to Nextflow |
| `--max-cpus` | `16` | CPU cap passed to Nextflow |
| `--nf-args` | _(none)_ | Extra pipeline flags, e.g. `"--orthofinder_v2"` |
| `--force` | off | Overwrite existing `run.sh` files |

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

1. Fix the accession in `benchmark/inputs/{genome_size}_{phylogeny}_{quality}.csv`
2. Regenerate the run directory: `./benchmark/run_benchmark.sh ... --force`
3. Re-run with `-resume` — Nextflow skips already-completed tasks:
   ```bash
   cd benchmark/results/bacteria_close_contiguous_n10 && ./run.sh -resume
   ```

### Nextflow session lock errors

If you see `Unable to acquire lock on session`, a previous run was interrupted
and left a stale lock. Clear it with:

```bash
RUN_ID=bacteria_close_contiguous_n10
# Check what's holding the lock (should be empty if nothing is running)
lsof benchmark/results/${RUN_ID}/cache/cache/*/db/LOCK

# Delete the stale cache
rm -rf benchmark/results/${RUN_ID}/cache

# Re-run fresh (no -resume)
cd benchmark/results/${RUN_ID} && ./run.sh
```

Each run directory has its own `cache/` folder so runs can never lock each other out.

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
