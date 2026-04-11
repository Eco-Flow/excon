#!/usr/bin/env python3
"""
collect_metrics.py — Parse EXCON Nextflow execution traces and produce
a benchmark comparison table suitable for a methods paper.

Scans benchmark/results/ for completed runs (directories containing
output/pipeline_info/execution_trace.tsv) and aggregates their metrics.

Usage:
    python3 collect_metrics.py \
        --results-dir benchmark/results \
        --output      benchmark/results/benchmark_metrics.tsv

Outputs:
    benchmark_metrics.tsv      — one row per run, all summary metrics
    benchmark_per_process.tsv  — one row per (run × process), detailed breakdown
"""

import argparse
import csv
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

_PER_SPECIES = re.compile(
    r'^(NCBIGENOMEDOWNLOAD|GUNZIP|AGAT_SPKEEPLONGESTISOFORM|GFFREAD|RENAME_FASTA)\s*\(.*\)$'
)
_CAFE_RUN_K = re.compile(r'^CAFE_RUN_K?\s*\(.*\)$', re.IGNORECASE)


def normalise_process(name: str) -> str:
    if _PER_SPECIES.match(name):
        return name.split('(')[0].strip()
    if _CAFE_RUN_K.match(name):
        return 'CAFE_RUN_K'
    return re.sub(r'\s*\(.*\)$', '', name).strip()


# ---------------------------------------------------------------------------
# Run ID parser — infer metadata from directory name
# ---------------------------------------------------------------------------

_RUN_ID_RE = re.compile(
    r'^(\w+)_(close|diverse)_(contiguous|fragmented)_n(\d+)$'
)

def parse_run_id(run_id: str) -> Optional[dict]:
    """Return {genome_size, phylogeny, quality, n_species} or None if unrecognised."""
    m = _RUN_ID_RE.match(run_id)
    if not m:
        return None
    return {
        'genome_size': m.group(1),
        'phylogeny':   m.group(2),
        'quality':     m.group(3),
        'n_species':   int(m.group(4)),
    }


# ---------------------------------------------------------------------------
# Trace file parser
# ---------------------------------------------------------------------------

def parse_trace(trace_path: Path) -> list:
    tasks = []
    with open(trace_path, newline='') as fh:
        reader = csv.DictReader(fh, delimiter='\t')
        for row in reader:
            if row.get('status', '') not in ('COMPLETED', 'CACHED'):
                continue
            tasks.append({
                'raw_name':    row.get('name', ''),
                'process':     normalise_process(row.get('name', '')),
                'status':      row.get('status', ''),
                'submit':      parse_datetime(row.get('submit', '')),
                'duration_s':  parse_duration(row.get('duration', '')),
                'realtime_s':  parse_duration(row.get('realtime', '')),
                'cpu_frac':    parse_cpu_pct(row.get('%cpu', '')),
                'peak_rss_b':  parse_size_bytes(row.get('peak_rss', '')),
                'peak_vmem_b': parse_size_bytes(row.get('peak_vmem', '')),
            })
    return tasks


# ---------------------------------------------------------------------------
# Per-run metric calculation
# ---------------------------------------------------------------------------

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
    'EGGNOG_DOWNLOAD',
    'EGGNOGMAPPER',
    'EGGNOG_TO_GO',
    'EGGNOG_TO_OG_GO',
    'OG_ANNOTATION_SUMMARY',
    'CAFE_GO_PREP',
    'CAFE_GO_RUN',
]


def compute_run_metrics(tasks: list, n_species: int, max_cpus: int) -> dict:
    if not tasks:
        return {}

    submits   = [t['submit'] for t in tasks if t['submit']]
    end_times = [
        t['submit'] + timedelta(seconds=t['duration_s'])
        for t in tasks if t['submit']
    ]
    total_wall_s = (max(end_times) - min(submits)).total_seconds() if submits else 0.0
    total_wall_h = total_wall_s / 3600

    cpu_seconds = sum(t['realtime_s'] * t['cpu_frac'] for t in tasks)
    cpu_hours   = cpu_seconds / 3600

    peak_ram_bytes = max((t['peak_rss_b'] for t in tasks), default=0)
    peak_ram_gb    = peak_ram_bytes / (1024 ** 3)

    allocated_cpu_hours = max_cpus * total_wall_h
    cpu_efficiency = cpu_hours / allocated_cpu_hours if allocated_cpu_hours > 0 else 0.0

    total_realtime = sum(t['realtime_s'] for t in tasks)
    mean_cpu_pct = (
        sum(t['cpu_frac'] * t['realtime_s'] for t in tasks) / total_realtime * 100
        if total_realtime > 0 else 0.0
    )

    genomes_per_hour = n_species / total_wall_h if total_wall_h > 0 else 0.0

    return {
        'total_wall_time_s':   round(total_wall_s, 1),
        'total_wall_time_min': round(total_wall_s / 60, 2),
        'total_cpu_hours':     round(cpu_hours, 4),
        'peak_ram_gb':         round(peak_ram_gb, 3),
        'mean_cpu_pct':        round(mean_cpu_pct, 1),
        'cpu_efficiency':      round(cpu_efficiency, 4),
        'genomes_per_hour':    round(genomes_per_hour, 4),
        'n_tasks':             len(tasks),
    }


def compute_per_process(tasks: list) -> dict:
    groups: dict = defaultdict(list)
    for t in tasks:
        groups[t['process']].append(t)

    result = {}
    for proc, ts in groups.items():
        realtime_total = sum(t['realtime_s'] for t in ts)
        submits   = [t['submit'] for t in ts if t['submit']]
        end_times = [
            t['submit'] + timedelta(seconds=t['duration_s'])
            for t in ts if t['submit']
        ]
        wall_s = (
            (max(end_times) - min(submits)).total_seconds()
            if submits else sum(t['duration_s'] for t in ts)
        )
        cpu_seconds = sum(t['realtime_s'] * t['cpu_frac'] for t in ts)
        peak_rss_gb = max(t['peak_rss_b'] for t in ts) / (1024 ** 3)
        mean_cpu_pct = (
            sum(t['cpu_frac'] * t['realtime_s'] for t in ts) / realtime_total * 100
            if realtime_total > 0 else 0.0
        )
        result[proc] = {
            'n_tasks':      len(ts),
            'wall_time_s':  round(wall_s, 1),
            'cpu_hours':    round(cpu_seconds / 3600, 4),
            'peak_ram_gb':  round(peak_rss_gb, 3),
            'mean_cpu_pct': round(mean_cpu_pct, 1),
        }
    return result


# ---------------------------------------------------------------------------
# Main
# ---------------------------------------------------------------------------

def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument('--results-dir', required=True,
                        help='benchmark/results directory containing run subdirs')
    parser.add_argument('--output', required=True,
                        help='Path for the summary TSV output')
    parser.add_argument('--max-cpus', type=int, default=16,
                        help='Allocated CPUs (for efficiency calculation)')
    args = parser.parse_args()

    results_dir  = Path(args.results_dir)
    summary_path = Path(args.output)
    per_proc_path = summary_path.with_name('benchmark_per_process.tsv')

    # Find the latest trace file per run directory.
    # benchmark.config writes timestamped filenames (execution_trace_YYYY-MM-DD_HH-mm-ss.tsv)
    # so each resume produces a new file; we take the most recent one per run.
    trace_files = []
    for run_dir in sorted(results_dir.iterdir()):
        if not run_dir.is_dir():
            continue
        pipeline_info = run_dir / 'output' / 'pipeline_info'
        if not pipeline_info.is_dir():
            continue
        candidates = sorted(pipeline_info.glob('execution_trace*.tsv'))
        if candidates:
            trace_files.append(candidates[-1])  # alphabetical = chronological for ISO timestamps

    if not trace_files:
        print("No completed runs found (no execution_trace*.tsv files).", file=sys.stderr)
        print(f"Expected pattern: {results_dir}/<run_id>/output/pipeline_info/execution_trace*.tsv",
              file=sys.stderr)
        sys.exit(0)

    summary_rows     = []
    per_process_rows = []

    for trace_file in trace_files:
        run_id = trace_file.parts[-4]   # results/<run_id>/output/pipeline_info/...
        meta   = parse_run_id(run_id)
        if meta is None:
            print(f"WARNING: Cannot parse run_id '{run_id}' — skipping", file=sys.stderr)
            continue

        tasks = parse_trace(trace_file)
        if not tasks:
            print(f"WARNING: Empty/unparseable trace for {run_id} — skipping", file=sys.stderr)
            continue

        metrics  = compute_run_metrics(tasks, meta['n_species'], args.max_cpus)
        per_proc = compute_per_process(tasks)

        row = {'run_id': run_id, **meta, **metrics}
        for stage in CAFE_STAGE_ORDER:
            p = per_proc.get(stage, {})
            row[f'wall_{stage.lower()}_s'] = p.get('wall_time_s', '')
        summary_rows.append(row)

        for proc, p in per_proc.items():
            per_process_rows.append({'run_id': run_id, **meta, 'process': proc, **p})

    if summary_rows:
        fields = list(summary_rows[0].keys())
        with open(summary_path, 'w', newline='') as fh:
            w = csv.DictWriter(fh, fieldnames=fields, delimiter='\t')
            w.writeheader()
            w.writerows(summary_rows)
        print(f"Summary written: {summary_path}  ({len(summary_rows)} runs)")

    if per_process_rows:
        fields2 = list(per_process_rows[0].keys())
        with open(per_proc_path, 'w', newline='') as fh:
            w = csv.DictWriter(fh, fieldnames=fields2, delimiter='\t')
            w.writeheader()
            w.writerows(per_process_rows)
        print(f"Per-process breakdown: {per_proc_path}")

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
