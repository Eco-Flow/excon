#!/usr/bin/env python3
"""
Compare N0.tsv and species tree outputs from two OrthoFinder runs.
Usage:
    compare_orthofinder_runs.py \\
        --run1 /path/to/run1/N0.tsv \\
        --run2 /path/to/run2/N0.tsv \\
        --tree1 /path/to/run1/SpeciesTree_rooted_node_labels.txt \\
        --tree2 /path/to/run2/SpeciesTree_rooted_node_labels.txt \\
        [--out comparison_report.txt]
"""

import argparse
import sys
import re
from collections import defaultdict


# ---------------------------------------------------------------------------
# Tree utilities (no dependencies beyond stdlib)
# ---------------------------------------------------------------------------

def parse_branch_lengths(newick):
    """Return sorted list of all branch lengths in a Newick string."""
    return sorted(float(x) for x in re.findall(r':([0-9]+\.?[0-9]*(?:[eE][+-]?[0-9]+)?)', newick))


def tip_to_root_distances(newick):
    """
    Return dict of {tip_name: root_to_tip_distance} by summing branch
    lengths along the path from root to each tip. Uses a simple recursive
    Newick parser — no ete3 / ape needed.
    """
    # Strip trailing semicolon and whitespace
    s = newick.strip().rstrip(';')

    def parse(s, depth=0.0):
        """Returns list of (tip_name, total_depth) pairs."""
        s = s.strip()
        if not s.startswith('('):
            # Leaf node: "name:length" or just "name"
            if ':' in s:
                name, length = s.rsplit(':', 1)
                return [(name.strip(), depth + float(length))]
            else:
                return [(s.strip(), depth)]

        # Find matching closing paren
        level = 0
        split_at = []
        for i, c in enumerate(s):
            if c == '(':
                level += 1
            elif c == ')':
                level -= 1
                if level == 0:
                    close = i
                    break
            elif c == ',' and level == 1:
                split_at.append(i)

        inner = s[1:close]
        after = s[close + 1:]  # e.g. "N4:79.8" or ":79.8" or ""

        # Branch length for this internal node
        bl = 0.0
        if ':' in after:
            bl_str = after.split(':')[1].split(')')[0].split(',')[0]
            try:
                bl = float(bl_str)
            except ValueError:
                bl = 0.0

        # Split inner on top-level commas
        children = []
        prev = 0
        level = 0
        for i, c in enumerate(inner):
            if c == '(':
                level += 1
            elif c == ')':
                level -= 1
            elif c == ',' and level == 0:
                children.append(inner[prev:i])
                prev = i + 1
        children.append(inner[prev:])

        result = []
        for child in children:
            result.extend(parse(child.strip(), depth + bl))
        return result

    return dict(parse(s))


# ---------------------------------------------------------------------------
# N0.tsv utilities
# ---------------------------------------------------------------------------

def load_n0(path):
    """
    Load N0.tsv (or Orthogroups.tsv).
    Returns (id_col, species_list, hog_counts) where hog_counts is a dict:
        { hog_id: { species: count, ... } }
    """
    with open(path) as fh:
        header = fh.readline().rstrip('\n').split('\t')

    # Identify ID column and extra metadata columns to skip
    skip = {'HOG', 'OG', 'Gene Tree Parent Clade', 'Orthogroup', 'Desc'}
    id_col = 'HOG' if 'HOG' in header else 'Orthogroup'
    species = [c for c in header if c not in skip]

    hog_counts = {}
    with open(path) as fh:
        fh.readline()  # skip header
        for line in fh:
            parts = line.rstrip('\n').split('\t')
            row = dict(zip(header, parts))
            hog_id = row[id_col]
            counts = {}
            for sp in species:
                cell = row.get(sp, '').strip()
                if cell == '':
                    continue
                # Cell can be gene IDs (comma-separated) or a plain integer
                if cell.isdigit():
                    n = int(cell)
                else:
                    n = len([g for g in cell.split(',') if g.strip()])
                counts[sp] = n
            hog_counts[hog_id] = counts

    return id_col, species, hog_counts


def hog_stats(hog_counts, species):
    """
    For each HOG compute max, min, differential, n_species_present.
    Returns list of dicts.
    """
    rows = []
    for hid, counts in hog_counts.items():
        vals = [counts.get(sp, 0) for sp in species]
        n_present = sum(1 for v in vals if v > 0)
        if n_present < 2:
            continue  # would be filtered anyway
        rows.append({
            'hog': hid,
            'max': max(vals),
            'min': min(vals),
            'diff': max(vals) - min(vals),
            'mean': sum(vals) / len(vals),
            'n_present': n_present,
        })
    return rows


def summarise(stats, label, out):
    diffs = [r['diff'] for r in stats]
    maxes = [r['max'] for r in stats]
    n = len(diffs)
    diffs_s = sorted(diffs)

    def pct(lst, p):
        i = int(len(lst) * p / 100)
        return lst[min(i, len(lst) - 1)]

    print(f"\n{'='*60}", file=out)
    print(f"  {label}", file=out)
    print(f"{'='*60}", file=out)
    print(f"  HOGs (>=2 species):     {n:,}", file=out)
    print(f"  Max gene count (any):   {max(maxes):,}", file=out)
    print(f"  Size differential:", file=out)
    print(f"    min    = {min(diffs)}", file=out)
    print(f"    median = {pct(diffs_s, 50)}", file=out)
    print(f"    p90    = {pct(diffs_s, 90)}", file=out)
    print(f"    p95    = {pct(diffs_s, 95)}", file=out)
    print(f"    p99    = {pct(diffs_s, 99)}", file=out)
    print(f"    max    = {max(diffs)}", file=out)

    for thresh in [50, 100, 200]:
        n_above = sum(1 for d in diffs if d > thresh)
        print(f"  HOGs with diff > {thresh:3d}:  {n_above:,}  ({100*n_above/n:.1f}%)", file=out)

    print(f"\n  Top 10 highest-differential HOGs:", file=out)
    top10 = sorted(stats, key=lambda r: -r['diff'])[:10]
    print(f"  {'HOG':<25} {'min':>5} {'max':>5} {'diff':>5} {'n_sp':>5}", file=out)
    for r in top10:
        print(f"  {r['hog']:<25} {r['min']:>5} {r['max']:>5} {r['diff']:>5} {r['n_present']:>5}", file=out)


# ---------------------------------------------------------------------------
# Main
# ---------------------------------------------------------------------------

def main():
    ap = argparse.ArgumentParser(description='Compare two OrthoFinder N0.tsv outputs')
    ap.add_argument('--run1',  required=True, help='N0.tsv from run 1 (working)')
    ap.add_argument('--run2',  required=True, help='N0.tsv from run 2 (failing)')
    ap.add_argument('--tree1', help='SpeciesTree_rooted_node_labels.txt from run 1')
    ap.add_argument('--tree2', help='SpeciesTree_rooted_node_labels.txt from run 2')
    ap.add_argument('--out',   help='Output file (default: stdout)')
    args = ap.parse_args()

    out = open(args.out, 'w') if args.out else sys.stdout

    print("OrthoFinder run comparison", file=out)
    print("=" * 60, file=out)

    # --- Load N0 tables ---
    id1, sp1, hogs1 = load_n0(args.run1)
    id2, sp2, hogs2 = load_n0(args.run2)

    print(f"\nRun 1:  {args.run1}", file=out)
    print(f"  ID column : {id1}", file=out)
    print(f"  Species   : {len(sp1)}", file=out)
    print(f"  HOGs total: {len(hogs1):,}", file=out)

    print(f"\nRun 2:  {args.run2}", file=out)
    print(f"  ID column : {id2}", file=out)
    print(f"  Species   : {len(sp2)}", file=out)
    print(f"  HOGs total: {len(hogs2):,}", file=out)

    # Species set differences
    s1, s2 = set(sp1), set(sp2)
    if s1 != s2:
        print(f"\n  *** Species mismatch ***", file=out)
        print(f"  Only in run1: {sorted(s1 - s2)}", file=out)
        print(f"  Only in run2: {sorted(s2 - s1)}", file=out)
    else:
        print(f"\n  Species sets are identical.", file=out)

    common_sp = sorted(s1 & s2)

    # --- Per-run stats ---
    stats1 = hog_stats(hogs1, sp1)
    stats2 = hog_stats(hogs2, sp2)

    summarise(stats1, f"Run 1 - {args.run1}", out)
    summarise(stats2, f"Run 2 - {args.run2}", out)

    # --- Tree comparison ---
    if args.tree1 and args.tree2:
        print(f"\n{'='*60}", file=out)
        print(f"  Tree comparison", file=out)
        print(f"{'='*60}", file=out)

        with open(args.tree1) as fh:
            nwk1 = fh.read().strip()
        with open(args.tree2) as fh:
            nwk2 = fh.read().strip()

        bl1 = parse_branch_lengths(nwk1)
        bl2 = parse_branch_lengths(nwk2)

        def bl_summary(bls, label):
            print(f"\n  {label}", file=out)
            print(f"    branches : {len(bls)}", file=out)
            print(f"    min      : {min(bls):.6f}", file=out)
            print(f"    median   : {bls[len(bls)//2]:.6f}", file=out)
            print(f"    max      : {max(bls):.6f}", file=out)
            print(f"    sum      : {sum(bls):.4f}", file=out)
            print(f"    near-zero (< 0.001): {sum(1 for b in bls if b < 0.001)}", file=out)

        bl_summary(bl1, f"Run 1: {args.tree1}")
        bl_summary(bl2, f"Run 2: {args.tree2}")

        # Tip-to-root depths
        tips1 = tip_to_root_distances(nwk1)
        tips2 = tip_to_root_distances(nwk2)

        common_tips = sorted(set(tips1) & set(tips2))
        print(f"\n  Root-to-tip distances (common species):", file=out)
        print(f"  {'species':<40} {'run1':>10} {'run2':>10} {'delta':>10}", file=out)
        for sp in common_tips:
            d1 = tips1[sp]
            d2 = tips2[sp]
            print(f"  {sp:<40} {d1:>10.4f} {d2:>10.4f} {d2-d1:>+10.4f}", file=out)

        # Ultrametricity check
        if tips1:
            depths1 = list(tips1.values())
            range1 = max(depths1) - min(depths1)
            print(f"\n  Run 1 tip-depth range (0 = ultrametric): {range1:.6f}", file=out)
        if tips2:
            depths2 = list(tips2.values())
            range2 = max(depths2) - min(depths2)
            print(f"  Run 2 tip-depth range (0 = ultrametric): {range2:.6f}", file=out)

    if args.out:
        out.close()
        print(f"Report written to {args.out}")


if __name__ == '__main__':
    main()
