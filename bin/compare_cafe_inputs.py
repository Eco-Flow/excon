#!/usr/bin/env python3
"""
Compare the exact input files passed to CAFE5 in two runs.

Usage:
    compare_cafe_inputs.py \\
        --tree1 /path/run2/work/a5/5ac667.../pruned_tree \\
        --tree2 /path/run7/work/ed/c84b49.../SpeciesTree_rooted_ultra.txt \\
        --counts1 /path/run2/work/a5/5ac667.../hog_gene_counts.tsv \\
        --counts2 /path/run7/work/ed/c84b49.../hog_gene_counts.tsv \\
        [--out report.txt]
"""

import argparse
import sys
import re
import math


# ---------------------------------------------------------------------------
# Tree parsing
# ---------------------------------------------------------------------------

def parse_branch_lengths(newick):
    return [float(x) for x in re.findall(
        r':([0-9]+\.?[0-9]*(?:[eE][+-]?[0-9]+)?)', newick)]


def tip_root_distances(newick):
    """Return dict {tip_name: root_to_tip_distance}."""
    s = newick.strip().rstrip(';')

    def parse(s, depth=0.0):
        s = s.strip()
        if not s.startswith('('):
            if ':' in s:
                name, length = s.rsplit(':', 1)
                # strip trailing node labels like )N4
                length = re.split(r'[,)\s]', length)[0]
                try:
                    return [(name.strip(), depth + float(length))]
                except ValueError:
                    return [(name.strip(), depth)]
            return [(s.strip(), depth)]

        # find matching close paren
        level = 0
        close = -1
        for i, c in enumerate(s):
            if c == '(':
                level += 1
            elif c == ')':
                level -= 1
                if level == 0:
                    close = i
                    break

        inner = s[1:close]
        after = s[close + 1:]
        bl = 0.0
        m = re.match(r'[^:]*:([0-9]+\.?[0-9]*(?:[eE][+-]?[0-9]+)?)', after)
        if m:
            bl = float(m.group(1))

        children = []
        prev, level = 0, 0
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


def percentile(lst, p):
    if not lst:
        return float('nan')
    s = sorted(lst)
    i = int(len(s) * p / 100)
    return s[min(i, len(s) - 1)]


def tree_report(path, label, out):
    with open(path) as fh:
        nwk = fh.read().strip()

    bls = parse_branch_lengths(nwk)
    tips = tip_root_distances(nwk)
    depths = list(tips.values())

    print(f"\n{'='*60}", file=out)
    print(f"  Tree: {label}", file=out)
    print(f"  File: {path}", file=out)
    print(f"{'='*60}", file=out)
    print(f"  Branches          : {len(bls)}", file=out)
    print(f"  Tips              : {len(depths)}", file=out)
    print(f"  Branch lengths:", file=out)
    print(f"    min             : {min(bls):.6f}", file=out)
    print(f"    median          : {percentile(bls, 50):.6f}", file=out)
    print(f"    max             : {max(bls):.6f}", file=out)
    print(f"    near-zero(<0.001): {sum(1 for b in bls if b < 0.001)}", file=out)
    if depths:
        drange = max(depths) - min(depths)
        print(f"  Root-to-tip depths:", file=out)
        print(f"    min             : {min(depths):.6f}", file=out)
        print(f"    max             : {max(depths):.6f}", file=out)
        print(f"    range (0=ultrametric): {drange:.6f}", file=out)
        print(f"    % deviation     : {100*drange/max(depths):.2f}%", file=out)
    return tips, bls


# ---------------------------------------------------------------------------
# Gene count comparison
# ---------------------------------------------------------------------------

def load_counts(path):
    """Return (species_list, dict{hog: {sp: count}})"""
    with open(path) as fh:
        header = fh.readline().rstrip('\n').split('\t')

    skip = {'Desc', 'HOG', 'OG', 'Gene Tree Parent Clade', 'Orthogroup'}
    id_col = 'HOG' if 'HOG' in header else ('Orthogroup' if 'Orthogroup' in header else header[1])
    species = [c for c in header if c not in skip]

    hogs = {}
    with open(path) as fh:
        fh.readline()
        for line in fh:
            parts = line.rstrip('\n').split('\t')
            row = dict(zip(header, parts))
            hid = row.get(id_col, '').strip()
            if not hid:
                continue
            counts = {}
            for sp in species:
                cell = row.get(sp, '').strip()
                if cell == '' or cell == '0':
                    counts[sp] = 0
                elif cell.isdigit():
                    counts[sp] = int(cell)
                else:
                    counts[sp] = len([g for g in cell.split(',') if g.strip()])
            hogs[hid] = counts
    return species, hogs


def count_stats(hogs, species):
    rows = []
    for hid, counts in hogs.items():
        vals = [counts.get(sp, 0) for sp in species]
        n_present = sum(1 for v in vals if v > 0)
        if n_present < 2:
            continue
        rows.append({
            'hog': hid,
            'max': max(vals),
            'min': min(vals),
            'diff': max(vals) - min(vals),
            'n_zeros': sum(1 for v in vals if v == 0),
            'n_sp': len(species),
        })
    return rows


def counts_report(path, label, out):
    species, hogs = load_counts(path)
    stats = count_stats(hogs, species)
    diffs = [r['diff'] for r in stats]
    zeros = [r['n_zeros'] for r in stats]

    print(f"\n{'='*60}", file=out)
    print(f"  Gene counts: {label}", file=out)
    print(f"  File: {path}", file=out)
    print(f"{'='*60}", file=out)
    print(f"  Families (>=2 sp): {len(stats):,}", file=out)
    print(f"  Species           : {len(species)}", file=out)
    print(f"  Size differential:", file=out)
    print(f"    min    = {min(diffs)}", file=out)
    print(f"    p50    = {percentile(diffs, 50)}", file=out)
    print(f"    p90    = {percentile(diffs, 90)}", file=out)
    print(f"    p95    = {percentile(diffs, 95)}", file=out)
    print(f"    p99    = {percentile(diffs, 99)}", file=out)
    print(f"    max    = {max(diffs)}", file=out)
    for t in [50, 100]:
        n = sum(1 for d in diffs if d > t)
        print(f"  Families diff>{t}: {n:,}  ({100*n/len(diffs):.2f}%)", file=out)
    print(f"  Families with >=1 zero-count species:", file=out)
    n_zeros = sum(1 for z in zeros if z > 0)
    print(f"    {n_zeros:,}  ({100*n_zeros/len(stats):.1f}%)", file=out)
    print(f"  Mean zeros per family: {sum(zeros)/len(zeros):.2f}", file=out)

    print(f"\n  Top 10 highest-differential families:", file=out)
    print(f"  {'HOG':<25} {'min':>5} {'max':>5} {'diff':>5} {'zeros':>6}", file=out)
    for r in sorted(stats, key=lambda x: -x['diff'])[:10]:
        print(f"  {r['hog']:<25} {r['min']:>5} {r['max']:>5} {r['diff']:>5} {r['n_zeros']:>6}", file=out)

    return species, stats


# ---------------------------------------------------------------------------
# Tree vs counts species check
# ---------------------------------------------------------------------------

def check_species_match(tips, species, label, out):
    tip_names = set(tips.keys())
    sp_names = set(species)
    print(f"\n  Species in tree vs counts ({label}):", file=out)
    only_tree = tip_names - sp_names
    only_counts = sp_names - tip_names
    if only_tree:
        print(f"    Only in tree   : {sorted(only_tree)[:5]}{'...' if len(only_tree)>5 else ''}", file=out)
    if only_counts:
        print(f"    Only in counts : {sorted(only_counts)[:5]}{'...' if len(only_counts)>5 else ''}", file=out)
    if not only_tree and not only_counts:
        print(f"    Tree and counts species match exactly.", file=out)


# ---------------------------------------------------------------------------
# Main
# ---------------------------------------------------------------------------

def main():
    ap = argparse.ArgumentParser(description='Compare exact CAFE5 inputs from two runs')
    ap.add_argument('--tree1',   required=True, help='Tree file used in run 1 (working)')
    ap.add_argument('--tree2',   required=True, help='Tree file used in run 2 (failing)')
    ap.add_argument('--counts1', required=True, help='hog_gene_counts.tsv from run 1')
    ap.add_argument('--counts2', required=True, help='hog_gene_counts.tsv from run 2')
    ap.add_argument('--out',     help='Output file (default: stdout)')
    args = ap.parse_args()

    out = open(args.out, 'w') if args.out else sys.stdout

    print("CAFE5 input comparison: working run vs failing run", file=out)

    tips1, bls1 = tree_report(args.tree1,   "Run 1 (working)", out)
    tips2, bls2 = tree_report(args.tree2,   "Run 2 (failing)", out)

    sp1, stats1 = counts_report(args.counts1, "Run 1 (working)", out)
    sp2, stats2 = counts_report(args.counts2, "Run 2 (failing)", out)

    check_species_match(tips1, sp1, "Run 1", out)
    check_species_match(tips2, sp2, "Run 2", out)

    # Side-by-side tip depths for common species
    common = sorted(set(tips1) & set(tips2))
    if common:
        print(f"\n{'='*60}", file=out)
        print(f"  Root-to-tip depth comparison (common tip names)", file=out)
        print(f"{'='*60}", file=out)
        print(f"  {'species':<42} {'run1':>10} {'run2':>10} {'ratio':>8}", file=out)
        for sp in common:
            d1, d2 = tips1[sp], tips2[sp]
            ratio = d2 / d1 if d1 > 0 else float('inf')
            print(f"  {sp:<42} {d1:>10.4f} {d2:>10.4f} {ratio:>8.3f}", file=out)

    if args.out:
        out.close()
        print(f"Report written to {args.out}")


if __name__ == '__main__':
    main()
