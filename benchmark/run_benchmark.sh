#!/usr/bin/env bash
# =============================================================================
# EXCON Benchmark Runner
# Runs CAFE-only (no GO, no stats) across factor combinations:
#   Genome size:  bacteria | insect | mammal
#   Phylogeny:    close | diverse
#   Quality:      contiguous | fragmented
#   Dataset size: 10 | 50 | 100 (--dataset-sizes)
#
# Usage:
#   ./run_benchmark.sh [options]
#
# Options:
#   --genome-sizes    Comma-separated subset, e.g. bacteria,insect  [default: all]
#   --phylogenies     Comma-separated subset, e.g. close            [default: all]
#   --qualities       Comma-separated subset, e.g. contiguous       [default: all]
#   --dataset-sizes   Comma-separated, e.g. 10,50                   [default: 10,50,100]
#   --profile         Nextflow profile, e.g. singularity            [default: singularity]
#   --custom-config   Path to your HPC config file (SGE/SLURM etc) [default: none]
#   --max-memory      Memory cap per job                            [default: 128.GB]
#   --max-cpus        CPU cap per job                               [default: 16]
#   --no-resume       Disable -resume (forces fresh reruns)
#   --dry-run         Print commands without running them
# =============================================================================

set -euo pipefail

SCRIPT_DIR="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"
EXCON_DIR="$(dirname "$SCRIPT_DIR")"
INPUTS_DIR="${SCRIPT_DIR}/inputs"
RESULTS_DIR="${SCRIPT_DIR}/results"

# --- Defaults ---
GENOME_SIZES_ARR=(bacteria insect mammal)
DATASET_SIZES_ARR=(10 50 100)
QUALITIES_ARR=(contiguous fragmented)
PHYLOGENIES_ARR=(close diverse)
PROFILE="singularity"
CUSTOM_CONFIG=""
MAX_MEMORY="128.GB"
MAX_CPUS=16
RESUME=true
DRY_RUN=false
NXF_VER="${NXF_VER:-25.10.0}"

# --- Argument parsing ---
while [[ $# -gt 0 ]]; do
    case $1 in
        --genome-sizes)  IFS=',' read -ra GENOME_SIZES_ARR   <<< "$2"; shift 2 ;;
        --dataset-sizes) IFS=',' read -ra DATASET_SIZES_ARR  <<< "$2"; shift 2 ;;
        --qualities)     IFS=',' read -ra QUALITIES_ARR       <<< "$2"; shift 2 ;;
        --phylogenies)   IFS=',' read -ra PHYLOGENIES_ARR     <<< "$2"; shift 2 ;;
        --profile)       PROFILE="$2";      shift 2 ;;
        --custom-config) CUSTOM_CONFIG="$2"; shift 2 ;;
        --max-memory)    MAX_MEMORY="$2";   shift 2 ;;
        --max-cpus)      MAX_CPUS="$2";     shift 2 ;;
        --no-resume)     RESUME=false;      shift ;;
        --dry-run)       DRY_RUN=true;      shift ;;
        --help|-h)
            sed -n '/^# Usage/,/^# ===/{ /^# ===/d; s/^# \?//; p }' "$0"
            exit 0 ;;
        *) echo "Unknown argument: $1" >&2; exit 1 ;;
    esac
done

# --- Setup ---
mkdir -p "$RESULTS_DIR"
LOG="${RESULTS_DIR}/run_log.tsv"

if [[ ! -f "$LOG" ]]; then
    printf 'run_id\tgenome_size\tphylogeny\tquality\tn_species\tstatus\tstart_time\tend_time\ttrace_file\n' > "$LOG"
fi

RESUME_FLAG=""
[[ "$RESUME" == true ]] && RESUME_FLAG="-resume"

# --- Helper: parse "1m 18s" or "24h" or "892ms" → integer seconds ---
# (used only for pre-run estimation; actual timing comes from trace)

log_run() {
    local run_id="$1" genome_size="$2" phylogeny="$3" quality="$4"
    local n_species="$5" status="$6" start="$7" end="$8" trace="$9"
    printf '%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\n' \
        "$run_id" "$genome_size" "$phylogeny" "$quality" \
        "$n_species" "$status" "$start" "$end" "$trace" >> "$LOG"
}

# =============================================================================
# Main loop
# =============================================================================

total_runs=0
planned_runs=0

for genome_size in "${GENOME_SIZES_ARR[@]}"; do
    for phylogeny in "${PHYLOGENIES_ARR[@]}"; do
        for quality in "${QUALITIES_ARR[@]}"; do

            full_csv="${INPUTS_DIR}/${genome_size}_${phylogeny}_${quality}.csv"

            if [[ ! -f "$full_csv" ]]; then
                echo "WARNING: Input CSV not found: ${full_csv} — skipping"
                continue
            fi

            # Count non-empty lines (each is a species)
            total_species=$(grep -c '[^[:space:]]' "$full_csv" || true)

            for n_species in "${DATASET_SIZES_ARR[@]}"; do

                run_id="${genome_size}_${phylogeny}_${quality}_n${n_species}"
                outdir="${RESULTS_DIR}/${run_id}"
                workdir="${RESULTS_DIR}/${run_id}_work"
                cachedir="${RESULTS_DIR}/${run_id}_cache"
                subset_csv="${RESULTS_DIR}/${run_id}_input.csv"

                ((planned_runs++)) || true

                if [[ $n_species -gt $total_species ]]; then
                    echo "SKIP [${run_id}]: need ${n_species} species but only ${total_species} in CSV"
                    log_run "$run_id" "$genome_size" "$phylogeny" "$quality" \
                            "$n_species" "SKIPPED_INSUFFICIENT_SPECIES" "" "" ""
                    continue
                fi

                # Check if already successfully completed
                if grep -q "^${run_id}	.*	COMPLETED	" "$LOG" 2>/dev/null; then
                    echo "SKIP [${run_id}]: already logged as COMPLETED (use --no-resume to force rerun)"
                    continue
                fi

                echo ""
                echo "════════════════════════════════════════"
                echo " Run: ${run_id}  (${n_species}/${total_species} species)"
                echo "════════════════════════════════════════"

                # Subset the CSV to first n_species rows
                head -n "$n_species" "$full_csv" > "$subset_csv"

                # ---- Nextflow command ----
                NF_CMD=(
                    env NXF_VER="${NXF_VER}"
                    NXF_CACHE_DIR="${cachedir}"
                    nextflow run "${EXCON_DIR}/main.nf"
                    -profile "${PROFILE}"
                    --input "${subset_csv}"
                    --outdir "${outdir}"
                    --max_memory "${MAX_MEMORY}"
                    --max_cpus "${MAX_CPUS}"
                    -work-dir "${workdir}"
                    -with-trace "${outdir}/pipeline_info/execution_trace.tsv"
                    -with-timeline "${outdir}/pipeline_info/execution_timeline.html"
                    -with-report "${outdir}/pipeline_info/execution_report.html"
                )
                # Pass custom config (HPC scheduler, singularity cache, queue names, etc.)
                [[ -n "$CUSTOM_CONFIG" ]] && NF_CMD+=(--custom_config "${CUSTOM_CONFIG}")
                # Append resume flag if set
                [[ -n "$RESUME_FLAG" ]] && NF_CMD+=("$RESUME_FLAG")

                if [[ "$DRY_RUN" == true ]]; then
                    echo "  [DRY RUN] ${NF_CMD[*]}"
                    ((total_runs++)) || true
                    continue
                fi

                start_time=$(date -Iseconds)
                exit_status=0

                "${NF_CMD[@]}" 2>&1 | tee "${RESULTS_DIR}/${run_id}.log" || exit_status=$?

                end_time=$(date -Iseconds)

                # Find trace file
                trace_file="${outdir}/pipeline_info/execution_trace.tsv"
                if [[ ! -f "$trace_file" ]]; then
                    trace_file=$(ls "${outdir}/pipeline_info/execution_trace"*.txt 2>/dev/null | head -1 || echo "")
                fi

                if [[ $exit_status -eq 0 ]]; then
                    run_status="COMPLETED"
                    echo "OK  [${run_id}] completed successfully"
                else
                    run_status="FAILED"
                    echo "ERR [${run_id}] failed with exit ${exit_status} — check ${RESULTS_DIR}/${run_id}.log"
                fi

                log_run "$run_id" "$genome_size" "$phylogeny" "$quality" \
                        "$n_species" "$run_status" "$start_time" "$end_time" "$trace_file"

                ((total_runs++)) || true
            done
        done
    done
done

echo ""
echo "════════════════════════════════════════"
echo " All runs finished (${total_runs} executed)"
echo "════════════════════════════════════════"

if [[ "$DRY_RUN" == false && $total_runs -gt 0 ]]; then
    echo " Collecting metrics..."
    python3 "${SCRIPT_DIR}/collect_metrics.py" \
        --results-dir "${RESULTS_DIR}" \
        --log "${LOG}" \
        --output "${RESULTS_DIR}/benchmark_metrics.tsv"
    echo " Report written to: ${RESULTS_DIR}/benchmark_metrics.tsv"
    echo " Per-process breakdown: ${RESULTS_DIR}/benchmark_per_process.tsv"
fi
