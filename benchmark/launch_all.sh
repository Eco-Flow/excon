#!/usr/bin/env bash
# Launch every not-yet-started run.sh under a benchmark results directory,
# throttled to --max-concurrent Nextflow sessions running at once. run_benchmark.sh
# only writes run.sh files (see its own header comment) — this is what actually
# starts them, without submitting one Nextflow session (and its own burst of SGE
# child jobs, capped separately by the executor's queueSize in your
# --custom-config) per run_benchmark.sh-generated directory all at the same time.
#
# Each run.sh is launched in the foreground of its own backgrounded subshell
# (deliberately not run.sh's own -bg, which double-forks and detaches — this
# script needs to `wait` on it to know when a slot frees up), with output
# redirected to launch.log inside that run's own directory.
#
# Safe to re-run: a directory is skipped only if its run actually completed
# successfully, or a Nextflow process is currently running in it — anything
# else (failed, interrupted, never checked its own exit code) is relaunched.
# Merely checking whether output/ or .nextflow.log exist isn't enough: a run
# that failed early (e.g. a bad NCBI accession, an OOM'd task) still creates
# both, and would otherwise be silently treated as done forever.
#
# Usage:
#   ./benchmark/launch_all.sh [--max-concurrent N] [--results-dir DIR]
#
# Run this under nohup/screen/tmux — it stays alive for as long as it takes to
# work through every run, since a slot only frees up once that run.sh exits:
#   nohup ./benchmark/launch_all.sh --max-concurrent 3 > benchmark/launch_all.log 2>&1 &
#   disown

set -o pipefail
# Deliberately no `set -u`: empty-array expansion under nounset is inconsistent
# across bash versions (3.2, still the default on some systems, treats
# "${arr[@]}" on a zero-element array as an unbound-variable error; 4.4+
# doesn't), and this script relies on that expansion while pids is empty.

RESULTS_DIR="benchmark/results"
MAX_CONCURRENT=3

while [[ $# -gt 0 ]]; do
    case $1 in
        --max-concurrent) MAX_CONCURRENT="$2"; shift 2 ;;
        --results-dir)    RESULTS_DIR="$2"; shift 2 ;;
        -h|--help)
            sed -n '2,26p' "$0"
            exit 0
            ;;
        *) echo "Unknown option: $1" >&2; exit 1 ;;
    esac
done

if [[ ! -d "$RESULTS_DIR" ]]; then
    echo "ERROR: $RESULTS_DIR not found (run benchmark/run_benchmark.sh first)" >&2
    exit 1
fi

# True if a `nextflow run` process is currently running with its cwd inside
# $1 — matched by cwd (via /proc, Linux-only, fine on Myriad) rather than a
# PID this script itself tracked, so it also recognises runs started by hand
# or by an earlier, separate invocation of this script.
is_active() {
    local target
    target=$(cd "$1" && pwd) || return 1
    local pid
    for pid in $(pgrep -f 'nextflow run' 2>/dev/null); do
        [[ "$(readlink -f "/proc/$pid/cwd" 2>/dev/null)" == "$target" ]] && return 0
    done
    return 1
}

echo "$(date '+%H:%M:%S')  Starting — max $MAX_CONCURRENT concurrent Nextflow session(s), results dir: $RESULTS_DIR"

pids=()
launched=0
skipped=0

for d in "$RESULTS_DIR"/*/; do
    [[ -f "$d/run.sh" ]] || continue
    run_id=$(basename "$d")

    if [[ -f "$d/.nextflow.log" ]] && grep -q "Execution complete" "$d/.nextflow.log" 2>/dev/null; then
        echo "SKIP   [$run_id]: completed successfully"
        skipped=$((skipped + 1))
        continue
    fi

    if is_active "$d"; then
        echo "SKIP   [$run_id]: a Nextflow session is currently running in this directory"
        skipped=$((skipped + 1))
        continue
    fi

    if [[ -f "$d/.nextflow.log" ]]; then
        echo "$(date '+%H:%M:%S')  RETRY  [$run_id]: previous attempt neither completed nor is running — relaunching"
    fi

    # Throttle: poll for a free slot before launching the next one. Polling
    # rather than `wait -n` (bash 4.3+ only, and not guaranteed on every login
    # node) — the poll interval is negligible next to how long a run actually
    # takes.
    while [[ ${#pids[@]} -ge $MAX_CONCURRENT ]]; do
        sleep 5
        alive=()
        for pid in "${pids[@]}"; do
            kill -0 "$pid" 2>/dev/null && alive+=("$pid")
        done
        pids=("${alive[@]}")
    done

    echo "$(date '+%H:%M:%S')  LAUNCH [$run_id]  (${#pids[@]}/${MAX_CONCURRENT} slots were busy)"
    ( cd "$d" && ./run.sh > launch.log 2>&1 )   &
    pids+=("$!")
    launched=$((launched + 1))
    sleep 5   # stagger submissions slightly rather than firing in the same instant
done

echo "$(date '+%H:%M:%S')  All eligible runs launched (${launched} started, ${skipped} already running) — waiting for the last ${#pids[@]} to finish..."
wait
echo "$(date '+%H:%M:%S')  Done."
