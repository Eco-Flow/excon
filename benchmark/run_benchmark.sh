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
#   --parallel        Launch all runs simultaneously (recommended on HPC)
#   --no-resume       Disable -resume (forces fresh reruns)
#   --dry-run         Print commands without running them
#   --nf-args         Extra Nextflow/pipeline args (quoted), e.g. "--orthofinder_v2"
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
PARALLEL=false
RESUME=true
DRY_RUN=false
NF_EXTRA_ARGS=""
NXF_VER="${NXF_VER:-25.10.0}"

# --- Argument parsing ---
while [[ $# -gt 0 ]]; do
    case $1 in
        --genome-sizes)  IFS=',' read -ra GENOME_SIZES_ARR   <<< "$2"; shift 2 ;;
        --dataset-sizes) IFS=',' read -ra DATASET_SIZES_ARR  <<< "$2"; shift 2 ;;
        --qualities)     IFS=',' read -ra QUALITIES_ARR       <<< "$2"; shift 2 ;;
        --phylogenies)   IFS=',' read -ra PHYLOGENIES_ARR     <<< "$2"; shift 2 ;;
        --profile)       PROFILE="$2";       shift 2 ;;
        --custom-config) CUSTOM_CONFIG="$2"; shift 2 ;;
        --max-memory)    MAX_MEMORY="$2";    shift 2 ;;
        --max-cpus)      MAX_CPUS="$2";      shift 2 ;;
        --parallel)      PARALLEL=true;      shift ;;
        --no-resume)     RESUME=false;       shift ;;
        --dry-run)       DRY_RUN=true;       shift ;;
        --nf-args)       NF_EXTRA_ARGS="$2"; shift 2 ;;
        --help|-h)
            sed -n '/^# Usage/,/^# ===/{ /^# ===/d; s/^# \?//; p }' "$0"
            exit 0 ;;
        *) echo "Unknown argument: $1" >&2; exit 1 ;;
    esac
done

# --- Setup ---
mkdir -p "$RESULTS_DIR"
LOG="${RESULTS_DIR}/run_log.tsv"
LOG_LOCK="${RESULTS_DIR}/.run_log.lock"

if [[ ! -f "$LOG" ]]; then
    printf 'run_id\tgenome_size\tphylogeny\tquality\tn_species\tstatus\tstart_time\tend_time\ttrace_file\n' > "$LOG"
fi

RESUME_FLAG=""
[[ "$RESUME" == true ]] && RESUME_FLAG="-resume"

# --- Log helper (flock-protected so parallel writes don't corrupt the TSV) ---
log_run() {
    local run_id="$1" genome_size="$2" phylogeny="$3" quality="$4"
    local n_species="$5" status="$6" start="$7" end="$8" trace="$9"
    (
        flock 200
        printf '%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\n' \
            "$run_id" "$genome_size" "$phylogeny" "$quality" \
            "$n_species" "$status" "$start" "$end" "$trace" >> "$LOG"
    ) 200>"$LOG_LOCK"
}

# =============================================================================
# execute_run — runs one Nextflow benchmark and logs the result.
# Called directly (sequential) or as a background job (parallel).
# =============================================================================
execute_run() {
    local run_id="$1" genome_size="$2" phylogeny="$3" quality="$4"
    local n_species="$5" subset_csv="$6" outdir="$7" workdir="$8" cachedir="$9"

    echo ""
    echo "════════════════════════════════════════"
    echo " Run: ${run_id}  (${n_species} species)"
    echo "════════════════════════════════════════"

    local NF_CMD=(
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
    [[ -n "$CUSTOM_CONFIG" ]]  && NF_CMD+=(--custom_config "${CUSTOM_CONFIG}")
    [[ -n "$NF_EXTRA_ARGS" ]]  && read -ra _extra <<< "$NF_EXTRA_ARGS" && NF_CMD+=("${_extra[@]}")
    [[ -n "$RESUME_FLAG" ]]    && NF_CMD+=("$RESUME_FLAG")

    local start_time exit_status=0
    start_time=$(date -Iseconds)

    "${NF_CMD[@]}" 2>&1 | tee "${RESULTS_DIR}/${run_id}.log" || exit_status=$?

    local end_time run_status trace_file
    end_time=$(date -Iseconds)

    trace_file="${outdir}/pipeline_info/execution_trace.tsv"
    [[ ! -f "$trace_file" ]] && \
        trace_file=$(ls "${outdir}/pipeline_info/execution_trace"*.txt 2>/dev/null | head -1 || echo "")

    if [[ $exit_status -eq 0 ]]; then
        run_status="COMPLETED"
        echo "OK  [${run_id}] completed"
    else
        run_status="FAILED"
        echo "ERR [${run_id}] failed — check ${RESULTS_DIR}/${run_id}.log"
    fi

    log_run "$run_id" "$genome_size" "$phylogeny" "$quality" \
            "$n_species" "$run_status" "$start_time" "$end_time" "$trace_file"
}

# =============================================================================
# Main loop — build the list of runs, then execute sequentially or in parallel
# =============================================================================

declare -a PIDS=()
total_runs=0

for genome_size in "${GENOME_SIZES_ARR[@]}"; do
    for phylogeny in "${PHYLOGENIES_ARR[@]}"; do
        for quality in "${QUALITIES_ARR[@]}"; do

            full_csv="${INPUTS_DIR}/${genome_size}_${phylogeny}_${quality}.csv"

            if [[ ! -f "$full_csv" ]]; then
                echo "WARNING: Input CSV not found: ${full_csv} — skipping"
                continue
            fi

            total_species=$(grep -c '[^[:space:]]' "$full_csv" || true)

            for n_species in "${DATASET_SIZES_ARR[@]}"; do

                run_id="${genome_size}_${phylogeny}_${quality}_n${n_species}"
                outdir="${RESULTS_DIR}/${run_id}"
                workdir="${RESULTS_DIR}/${run_id}_work"
                cachedir="${RESULTS_DIR}/${run_id}_cache"
                subset_csv="${RESULTS_DIR}/${run_id}_input.csv"

                if [[ $n_species -gt $total_species ]]; then
                    echo "SKIP [${run_id}]: need ${n_species} species but only ${total_species} in CSV"
                    log_run "$run_id" "$genome_size" "$phylogeny" "$quality" \
                            "$n_species" "SKIPPED_INSUFFICIENT_SPECIES" "" "" ""
                    continue
                fi

                if grep -q "^${run_id}	.*	COMPLETED	" "$LOG" 2>/dev/null; then
                    echo "SKIP [${run_id}]: already COMPLETED"
                    continue
                fi

                # Remove pipeline_info from any previous failed run — Nextflow refuses
                # to start if execution_trace.tsv / execution_report.html already exist
                rm -rf "${outdir}/pipeline_info"

                # Prepare subset CSV (safe to do before backgrounding)
                head -n "$n_species" "$full_csv" > "$subset_csv"

                if [[ "$DRY_RUN" == true ]]; then
                    echo "  [DRY RUN] ${run_id}"
                    ((total_runs++)) || true
                    continue
                fi

                ((total_runs++)) || true

                if [[ "$PARALLEL" == true ]]; then
                    # Launch in background; stdout/stderr goes to per-run log
                    execute_run \
                        "$run_id" "$genome_size" "$phylogeny" "$quality" \
                        "$n_species" "$subset_csv" "$outdir" "$workdir" "$cachedir" \
                        >> "${RESULTS_DIR}/${run_id}.log" 2>&1 &
                    PIDS+=($!)
                    echo "LAUNCHED [${run_id}] — PID $!"
                else
                    execute_run \
                        "$run_id" "$genome_size" "$phylogeny" "$quality" \
                        "$n_species" "$subset_csv" "$outdir" "$workdir" "$cachedir"
                fi
            done
        done
    done
done

# --- Wait for all background jobs (parallel mode) ---
if [[ "$PARALLEL" == true && ${#PIDS[@]} -gt 0 ]]; then
    echo ""
    echo "Waiting for ${#PIDS[@]} parallel runs to complete..."
    for pid in "${PIDS[@]}"; do
        wait "$pid" || true   # errors are logged per-run; don't abort here
    done
fi

echo ""
echo "════════════════════════════════════════"
echo " All runs finished (${total_runs} launched)"
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
