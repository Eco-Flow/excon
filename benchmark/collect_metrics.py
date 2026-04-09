#!/usr/bin/env python3
"""
collect_metrics.py — Parse EXCON Nextflow execution traces and produce
a benchmark comparison table suitable for a methods paper.

Usage:
    python3 collect_metrics.py \
        --results-dir benchmark/results \
        --log         benchmark/results/run_log.tsv \
        --output      benchmark/results/benchmark_metrics.tsv

Outputs:
    benchmark_metrics.tsv      — one row per run, all summary metrics
    benchmark_per_process.tsv  — one row per (run × process), detailed breakdown
"""

import argparse
import csv
import os
import re
import sys
from collections import defaultdict
from datetime import datetime, timedelta
from pathlib import Path
from typing import Optional

# ---------------------------------------------------------------------------
# Time / size parsers — handle Nextflow's human-readable trace format
# ---------------------------------------------------------------------------

_TIME_UNITS = {
    'd':  86400,
    'h':  3600,
    'm':  60,
    's':  1,
    'ms': 0.001,
}

def parse_duration(s: str) -> float:
    """Return duration in seconds from strings like '3m 7s', '892ms', '1h 5m', '45s'."""
    if not s or s.strip() in ('-', ''):
        return 0.0
    s = s.strip()
    total = 0.0
    # Match tokens like "3m", "7s", "892ms", "1h", "1d"
    for value, unit in re.findall(r'(\d+(?:\.\d+)?)\s*(d|h|ms|m|s)', s):
        total += float(value) * _TIME_UNITS[unit]
    return total


def parse_size_bytes(s: str) -> float:
    """Return size in bytes from strings like '1.7 GB', '155.8 MB', '7.1 KB', '500 B'."""
    if not s or s.strip() in ('-', ''):
        return 0.0
    s = s.strip()
    m = re.match(r'([\d.]+)\s*([KMGTP]?B)', s, re.IGNORECASE)
    if not m:
        return 0.0
    value = float(m.group(1))
    unit = m.group(2).upper()
    multipliers = {'B': 1, 'KB': 1024, 'MB': 1024**2, 'GB': 1024**3, 'TB': 1024**4}
    return value * multipliers.get(unit, 1)


def parse_cpu_pct(s: str) -> float:
    """Return CPU as a fraction of one core from '12.3%' → 0.123, '572.8%' → 5.728."""
    if not s or s.strip() in ('-', ''):
        return 0.0
    return float(s.strip('%')) / 100.0


def parse_datetime(s: str) -> Optional[datetime]:
    """Parse '2026-04-08 21:32:49.530' → datetime."""
    if not s or s.strip() in ('-', ''):
        return None
    for fmt in ('%Y-%m-%d %H:%M:%S.%f', '%Y-%m-%d %H:%M:%S'):
        try:
            return datetime.strptime(s.strip(), fmt)
        except ValueError:
            pass
    return None


# ---------------------------------------------------------------------------
# Process-name normaliser — collapse per-species / per-k duplicates
# ---------------------------------------------------------------------------

# Processes that run once per species (strip the parenthetical suffix)
_PER_SPECIES = re.compile(
    r'^(NCBIGENOMEDOWNLOAD|GUNZIP|AGAT_SPKEEPLONGESTISOFORM|GFFREAD|RENAME_FASTA)\s*\(.*\)$'
)
# CAFE_RUN_K runs multiple times with different k values
_CAFE_RUN_K  = re.compile(r'^CAFE_RUN_K?\s*\(.*\)$', re.IGNORECASE)


def normalise_process(name: str) -> str:
    """Return a canonical process name for grouping."""
    if _PER_SPECIES.match(name):
        return name.split('(')[0].strip()
    if _CAFE_RUN_K.match(name):
        return 'CAFE_RUN_K'
    # Strip any trailing parenthetical for other processes
    return re.sub(r'\s*\(.*\)$', '', name).strip()


# ---------------------------------------------------------------------------
# Trace file parser
# ---------------------------------------------------------------------------

def parse_trace(trace_path: str) -> list[dict]:
    """Parse a Nextflow execution_trace TSV and return a list of task dicts."""
    tasks = []
    with open(trace_path, newline='') as fh:
        reader = csv.DictReader(fh, delimiter='\t')
        for row in reader:
            if row.get('status', '') not in ('COMPLETED', 'CACHED'):
                continue
            tasks.append({
                'raw_name':   row.get('name', ''),
                'process':    normalise_process(row.get('name', '')),
                'status':     row.get('status', ''),
                'submit':     parse_datetime(row.get('submit', '')),
                'duration_s': parse_duration(row.get('duration', '')),
                'realtime_s': parse_duration(row.get('realtime', '')),
                'cpu_frac':   parse_cpu_pct(row.get('%cpu', '')),
                'peak_rss_b': parse_size_bytes(row.get('peak_rss', '')),
                'peak_vmem_b':parse_size_bytes(row.get('peak_vmem', '')),
            })
    return tasks


# ---------------------------------------------------------------------------
# Per-run metric calculation
# ---------------------------------------------------------------------------

# Processes to include in the CAFE-pipeline breakdown (ordered for display)
CAFE_STAGE_ORDER = [
    'NCBIGENOMEDOWNLOAD',
    'GUNZIP',
    'AGAT_SPKEEPLONGESTISOFORM',
    'GFFREAD',
    'RENAME_FASTA',
    'ORTHOFINDER_CAFE',
    'ORTHOFINDER_V2_CAFE',
    'RESCALE_TREE',
    'CAFE_PREP',
    'CAFE_RUN_K',
    'CAFE_SELECT_K',
    'CAFE_RUN_BEST',
    'CAFE_MODEL_COMPARE',
    'CAFE_PLOT',
    'CAFE_RUN_LARGE',
]


def compute_run_metrics(tasks: list[dict], n_species: int, max_cpus: int) -> dict:
    """
    Compute summary metrics for one benchmark run.

    Returns a flat dict of scalar values.
    """
    if not tasks:
        return {}

    # --- Wall time = span from first submit to last submit + duration ---
    submits    = [t['submit'] for t in tasks if t['submit']]
    end_times  = [
        t['submit'] + timedelta(seconds=t['duration_s'])
        for t in tasks if t['submit']
    ]
    if not submits:
        total_wall_s = 0.0
    else:
        total_wall_s = (max(end_times) - min(submits)).total_seconds()

    total_wall_h = total_wall_s / 3600

    # --- CPU-seconds = sum(realtime_s × cpu_frac) across all tasks ---
    cpu_seconds = sum(t['realtime_s'] * t['cpu_frac'] for t in tasks)
    cpu_hours   = cpu_seconds / 3600

    # --- Peak RAM across whole run ---
    peak_ram_bytes = max((t['peak_rss_b'] for t in tasks), default=0)
    peak_ram_gb    = peak_ram_bytes / (1024 ** 3)

    # --- CPU efficiency = cpu_hours / (allocated_cpus × wall_hours) ---
    allocated_cpu_hours = max_cpus * total_wall_h
    cpu_efficiency = cpu_hours / allocated_cpu_hours if allocated_cpu_hours > 0 else 0.0

    # --- Mean %CPU across tasks (weighted by realtime) ---
    total_realtime = sum(t['realtime_s'] for t in tasks)
    if total_realtime > 0:
        mean_cpu_pct = sum(t['cpu_frac'] * t['realtime_s'] for t in tasks) / total_realtime * 100
    else:
        mean_cpu_pct = 0.0

    # --- Throughput ---
    genomes_per_hour = n_species / total_wall_h if total_wall_h > 0 else 0.0

    return {
        'total_wall_time_s':    round(total_wall_s, 1),
        'total_wall_time_min':  round(total_wall_s / 60, 2),
        'total_cpu_hours':      round(cpu_hours, 4),
        'peak_ram_gb':          round(peak_ram_gb, 3),
        'mean_cpu_pct':         round(mean_cpu_pct, 1),
        'cpu_efficiency':       round(cpu_efficiency, 4),
        'genomes_per_hour':     round(genomes_per_hour, 4),
        'n_tasks':              len(tasks),
    }


def compute_per_process(tasks: list[dict]) -> dict[str, dict]:
    """
    Aggregate tasks by process name.

    Returns {process_name: {wall_time_s, cpu_hours, peak_ram_gb, mean_cpu_pct, n_tasks}}.
    """
    groups: dict[str, list] = defaultdict(list)
    for t in tasks:
        groups[t['process']].append(t)

    result = {}
    for proc, ts in groups.items():
        realtime_total = sum(t['realtime_s'] for t in ts)
        # Wall time for a group = from first submit to last end
        submits   = [t['submit'] for t in ts if t['submit']]
        end_times = [
            t['submit'] + timedelta(seconds=t['duration_s'])
            for t in ts if t['submit']
        ]
        if submits:
            wall_s = (max(end_times) - min(submits)).total_seconds()
        else:
            wall_s = sum(t['duration_s'] for t in ts)

        cpu_seconds  = sum(t['realtime_s'] * t['cpu_frac'] for t in ts)
        peak_rss_gb  = max(t['peak_rss_b'] for t in ts) / (1024**3)
        mean_cpu_pct = (
            sum(t['cpu_frac'] * t['realtime_s'] for t in ts) / realtime_total * 100
            if realtime_total > 0 else 0.0
        )

        result[proc] = {
            'n_tasks':       len(ts),
            'wall_time_s':   round(wall_s, 1),
            'cpu_hours':     round(cpu_seconds / 3600, 4),
            'peak_ram_gb':   round(peak_rss_gb, 3),
            'mean_cpu_pct':  round(mean_cpu_pct, 1),
        }
    return result


# ---------------------------------------------------------------------------
# Main
# ---------------------------------------------------------------------------

def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument('--results-dir', required=True)
    parser.add_argument('--log',         required=True,
                        help='run_log.tsv produced by run_benchmark.sh')
    parser.add_argument('--output',      required=True,
                        help='Path for the summary TSV output')
    parser.add_argument('--max-cpus',    type=int, default=16,
                        help='Allocated CPUs (for efficiency calculation)')
    args = parser.parse_args()

    results_dir = Path(args.results_dir)
    summary_path    = Path(args.output)
    per_proc_path   = summary_path.with_name('benchmark_per_process.tsv')

    # --- Read run log ---
    runs = []
    with open(args.log, newline='') as fh:
        for row in csv.DictReader(fh, delimiter='\t'):
            if row.get('status') not in ('COMPLETED',):
                continue
            runs.append(row)

    if not runs:
        print("No COMPLETED runs found in log. Nothing to report.", file=sys.stderr)
        sys.exit(0)

    summary_rows    = []
    per_process_rows = []

    for run in runs:
        run_id    = run['run_id']
        n_species = int(run.get('n_species', 0))

        # Find trace file
        trace_file = run.get('trace_file', '').strip()
        if not trace_file or not os.path.exists(trace_file):
            # Fallback: search in outdir
            outdir = results_dir / run_id
            candidates = list(outdir.glob('pipeline_info/execution_trace*.ts*')) + \
                         list(outdir.glob('pipeline_info/execution_trace*.txt'))
            trace_file = str(candidates[0]) if candidates else ''

        if not trace_file or not os.path.exists(trace_file):
            print(f"WARNING: No trace found for {run_id} — skipping", file=sys.stderr)
            continue

        tasks = parse_trace(trace_file)
        if not tasks:
            print(f"WARNING: Empty/unparseable trace for {run_id} — skipping", file=sys.stderr)
            continue

        metrics = compute_run_metrics(tasks, n_species, args.max_cpus)
        per_proc = compute_per_process(tasks)

        # --- Summary row ---
        row = {
            'run_id':          run_id,
            'genome_size':     run.get('genome_size', ''),
            'phylogeny':       run.get('phylogeny', ''),
            'quality':         run.get('quality', ''),
            'n_species':       n_species,
            **metrics,
        }
        # Add per-stage wall times as extra columns (useful for paper tables)
        for stage in CAFE_STAGE_ORDER:
            p = per_proc.get(stage, {})
            row[f'wall_{stage.lower()}_s'] = p.get('wall_time_s', '')
        summary_rows.append(row)

        # --- Per-process rows ---
        for proc, p in per_proc.items():
            per_process_rows.append({
                'run_id':       run_id,
                'genome_size':  run.get('genome_size', ''),
                'phylogeny':    run.get('phylogeny', ''),
                'quality':      run.get('quality', ''),
                'n_species':    n_species,
                'process':      proc,
                **p,
            })

    # --- Write summary TSV ---
    if summary_rows:
        fields = list(summary_rows[0].keys())
        with open(summary_path, 'w', newline='') as fh:
            w = csv.DictWriter(fh, fieldnames=fields, delimiter='\t')
            w.writeheader()
            w.writerows(summary_rows)
        print(f"Summary written: {summary_path}  ({len(summary_rows)} runs)")

    # --- Write per-process TSV ---
    if per_process_rows:
        fields2 = list(per_process_rows[0].keys())
        with open(per_proc_path, 'w', newline='') as fh:
            w = csv.DictWriter(fh, fieldnames=fields2, delimiter='\t')
            w.writeheader()
            w.writerows(per_process_rows)
        print(f"Per-process breakdown: {per_proc_path}")

    # --- Print a quick console table ---
    if summary_rows:
        print()
        print(f"{'Run ID':<45} {'N':>4} {'Wall(min)':>10} {'CPU-h':>8} {'PeakRAM':>10} {'CPU-eff':>8} {'Gen/h':>8}")
        print('-' * 100)
        for r in summary_rows:
            print(
                f"{r['run_id']:<45}"
                f" {r['n_species']:>4}"
                f" {r.get('total_wall_time_min', ''):>10}"
                f" {r.get('total_cpu_hours', ''):>8}"
                f" {r.get('peak_ram_gb', ''):>9}GB"
                f" {r.get('cpu_efficiency', ''):>8}"
                f" {r.get('genomes_per_hour', ''):>8}"
            )


if __name__ == '__main__':
    main()
