#!/usr/bin/env python3

"""
Merge per-family CAFE5 output directories into one combined result directory.

CAFE_RUN_LARGE now runs each high-differential family on its own (see
cafe_run_large.nf) so it can fit whatever lambda that one family needs, instead
of forcing one lambda to explain every excluded family at once — which is what
made this module never converge. This script stitches those many single-family
runs back into one directory shaped like an ordinary CAFE5 result, so the
existing downstream consumers (cafeplotter via CAFE_PLOT_LARGE, CAFE_GO_PREP_LARGE,
CAFE_SIG_FAMILIES) can read it exactly as they would a normal run's output.

Every file but one is simple concatenation: each per-family run wrote a header
plus its own one data row. *_clade_results.txt is the exception — CAFE5 reports
it as a total across every family in the run, so a single-family run's file is
really just that family's own 0/1 contribution per node, and the per-family
files have to be summed rather than concatenated.

Pure standard library, no BioPython/pandas, to match the rest of bin/.
"""

import argparse
import csv
import glob
import os
import re
import sys


def find_one(dirpath, suffix):
    hits = sorted(glob.glob(os.path.join(dirpath, "*" + suffix)))
    return hits[0] if hits else None


def merge_concat_table(family_dirs, suffix, out_path):
    """Concatenate '<header>\\n<one row per family>' tables: *_change.tab,
    *_branch_probabilities.tab, *_count.tab, *_family_results.txt."""
    header = None
    rows = []
    for d in family_dirs:
        f = find_one(d, suffix)
        if f is None:
            continue
        with open(f) as fh:
            lines = [line.rstrip("\n") for line in fh if line.strip()]
        if not lines:
            continue
        if header is None:
            header = lines[0]
        rows.extend(lines[1:])
    if header is None:
        return 0
    with open(out_path, "w") as out:
        out.write(header + "\n")
        for row in rows:
            out.write(row + "\n")
    return len(rows)


def merge_clade_results(family_dirs, suffix, out_path):
    """Sum *_clade_results.txt (taxon, increase, decrease) across families."""
    header = None
    order = []
    totals = {}
    for d in family_dirs:
        f = find_one(d, suffix)
        if f is None:
            continue
        with open(f) as fh:
            rows = list(csv.reader(fh, delimiter="\t"))
        if not rows:
            continue
        if header is None:
            header = rows[0]
        for row in rows[1:]:
            if not row:
                continue
            taxon, inc, dec = row[0], int(row[1]), int(row[2])
            if taxon not in totals:
                totals[taxon] = [0, 0]
                order.append(taxon)
            totals[taxon][0] += inc
            totals[taxon][1] += dec
    if header is None:
        return 0
    with open(out_path, "w", newline="") as out:
        writer = csv.writer(out, delimiter="\t")
        writer.writerow(header)
        for taxon in order:
            inc, dec = totals[taxon]
            writer.writerow([taxon, inc, dec])
    return len(order)


def merge_asr_tree(family_dirs, suffix, out_path):
    """Splice per-family NEXUS *_asr.tre files into one multi-tree file.

    Each per-family run's NEXUS file has exactly one TREE line (it only ever
    contained that one family). The preamble/translate block is identical
    across families (same fixed species tree in every run), so it — and the
    closing lines — are taken from the first file found, with every family's
    own TREE line inserted in between.
    """
    preamble, tail = None, None
    tree_lines = []
    for d in family_dirs:
        f = find_one(d, suffix)
        if f is None:
            continue
        with open(f) as fh:
            lines = fh.readlines()
        tree_idx = [i for i, line in enumerate(lines) if re.match(r"^\s*TREE\b", line, re.IGNORECASE)]
        if not tree_idx:
            continue
        i = tree_idx[0]
        if preamble is None:
            preamble = lines[:i]
            tail = lines[i + 1:]
        tree_lines.append(lines[i])
    if preamble is None:
        return 0
    with open(out_path, "w") as out:
        out.writelines(preamble)
        out.writelines(tree_lines)
        out.writelines(tail)
    return len(tree_lines)


def parse_base_results(dirpath):
    """Pull this family's own fitted lambda and -lnL out of its Base_results.txt.

    Deliberately not a "*_results.txt" glob: that pattern also matches
    *_family_results.txt and *_clade_results.txt, which end the same way.
    """
    f = None
    for fname in ("Base_results.txt", "Gamma_results.txt"):
        candidate = os.path.join(dirpath, fname)
        if os.path.isfile(candidate):
            f = candidate
            break
    if f is None:
        return None, None
    lam, score = None, None
    with open(f) as fh:
        for line in fh:
            m = re.match(r"Lambda:\s*([0-9.eE+-]+)", line)
            if m:
                lam = m.group(1)
            m = re.match(r"Model \w+ Final Likelihood \(-lnL\):\s*([0-9.eE+-]+|inf)", line)
            if m:
                score = m.group(1)
    return lam, score


def main():
    ap = argparse.ArgumentParser(description=__doc__)
    ap.add_argument("--outdir", required=True, help="Merged CAFE5-shaped output directory to create")
    ap.add_argument("--lambda-summary", default=None,
                    help="Optional TSV recording each family's own fitted lambda/-lnL")
    ap.add_argument("family_dirs", nargs="+", help="Per-family CAFE5 result directories")
    args = ap.parse_args()

    family_dirs = sorted(d for d in args.family_dirs if os.path.isdir(d))
    if not family_dirs:
        sys.exit("ERROR: no per-family result directories supplied — nothing converged")

    os.makedirs(args.outdir, exist_ok=True)

    n_asr   = merge_asr_tree(family_dirs, "_asr.tre", os.path.join(args.outdir, "Base_asr.tre"))
    n_bp    = merge_concat_table(family_dirs, "_branch_probabilities.tab", os.path.join(args.outdir, "Base_branch_probabilities.tab"))
    n_chg   = merge_concat_table(family_dirs, "_change.tab", os.path.join(args.outdir, "Base_change.tab"))
    n_cnt   = merge_concat_table(family_dirs, "_count.tab", os.path.join(args.outdir, "Base_count.tab"))
    n_fam   = merge_concat_table(family_dirs, "_family_results.txt", os.path.join(args.outdir, "Base_family_results.txt"))
    n_clade = merge_clade_results(family_dirs, "_clade_results.txt", os.path.join(args.outdir, "Base_clade_results.txt"))

    print("Merged %d per-family CAFE5 runs into %s" % (len(family_dirs), args.outdir))
    print("  asr trees: %d  branch_probabilities rows: %d  change rows: %d"
          % (n_asr, n_bp, n_chg))
    print("  count rows: %d  family_results rows: %d  clade taxa: %d"
          % (n_cnt, n_fam, n_clade))

    if n_asr == 0:
        sys.exit("ERROR: none of the per-family runs produced a usable *_asr.tre — nothing converged")

    if args.lambda_summary:
        with open(args.lambda_summary, "w") as out:
            out.write("HOG\tlambda\tneg_lnL\n")
            for d in family_dirs:
                hog = os.path.basename(d.rstrip("/")).replace("Out_cafe_large_", "", 1)
                lam, score = parse_base_results(d)
                out.write("%s\t%s\t%s\n" % (hog, lam or "NA", score or "NA"))
        print("Per-family lambda summary written to: %s" % args.lambda_summary)


if __name__ == "__main__":
    main()
