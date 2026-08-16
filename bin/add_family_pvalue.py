#!/usr/bin/env python3

"""
Join CAFE5's family-wide p-value onto the branch-level significant-changes
table, so both cutoffs requested for the manuscript can be applied: family-wide
p < --family-pcut (this family's rate differs from the background model) AND
branch p <= --branch-pcut (this specific branch drove the change).

Post-hoc only — reads outputs a finished run already published under
--outdir, joins by HOG, and writes two new files alongside them. Does not
touch or re-run any pipeline step, so it can be pointed at any past or future
--outdir without changes to main.nf.

Reads (relative to --outdir):
  cafe/significant_families/combined_changes_per_node.tsv
  cafe/model_comparison/best_cafe_results/{Base,Gamma}_family_results.txt
  cafe/large_families/Base_family_results.txt

Writes (relative to --outdir):
  cafe/significant_families/combined_changes_per_node.with_family_pvalue.tsv
    — every changed branch, plus family_pvalue and dual_significant columns
  cafe/significant_families/combined_significant_changes_per_node.dual_cutoff.tsv
    — only the rows where dual_significant is True

Pure standard library, to match the rest of bin/.
"""

import argparse
import csv
import os
import sys


def read_family_pvalues(path):
    """FamilyID -> pvalue from a CAFE5 *_family_results.txt (#FamilyID, pvalue, ...)."""
    pvals = {}
    if not path or not os.path.isfile(path):
        return pvals
    with open(path, newline="") as fh:
        reader = csv.reader(fh, delimiter="\t")
        next(reader, None)  # header
        for row in reader:
            if len(row) < 2:
                continue
            try:
                pvals[row[0]] = float(row[1])
            except ValueError:
                continue
    return pvals


def first_existing(*paths):
    for p in paths:
        if os.path.isfile(p):
            return p
    return None


def main():
    ap = argparse.ArgumentParser(
        description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter
    )
    ap.add_argument("--outdir", required=True, help="pipeline --outdir of the finished run")
    ap.add_argument("--family-pcut", type=float, default=0.01,
                     help="family-wide p-value cutoff (default 0.01)")
    ap.add_argument("--branch-pcut", type=float, default=0.05,
                     help="branch p-value cutoff (default 0.05)")
    args = ap.parse_args()

    changes_path = os.path.join(
        args.outdir, "cafe", "significant_families", "combined_changes_per_node.tsv"
    )
    if not os.path.isfile(changes_path):
        sys.exit(f"ERROR: {changes_path} not found — is --outdir a finished run?")

    main_dir = os.path.join(args.outdir, "cafe", "model_comparison", "best_cafe_results")
    large_dir = os.path.join(args.outdir, "cafe", "large_families")

    main_fam_file = first_existing(
        os.path.join(main_dir, "Base_family_results.txt"),
        os.path.join(main_dir, "Gamma_family_results.txt"),
    )
    large_fam_file = first_existing(os.path.join(large_dir, "Base_family_results.txt"))

    main_pvals = read_family_pvalues(main_fam_file)
    large_pvals = read_family_pvalues(large_fam_file)
    print(f"Loaded {len(main_pvals)} main-model family p-values from {main_fam_file}")
    print(f"Loaded {len(large_pvals)} large-family p-values from {large_fam_file}")

    with open(changes_path, newline="") as fh:
        reader = csv.reader(fh, delimiter="\t")
        header = next(reader)
        idx = {name: i for i, name in enumerate(header)}
        for required in ("HOG", "pvalue", "source"):
            if required not in idx:
                sys.exit(f"ERROR: {changes_path} has no '{required}' column")
        rows = list(reader)

    missing = 0
    out_rows = []
    for row in rows:
        hog = row[idx["HOG"]]
        source = row[idx["source"]]
        table = large_pvals if source == "large_family_own_lambda" else main_pvals
        fam_p = table.get(hog)
        if fam_p is None:
            missing += 1

        try:
            branch_p = float(row[idx["pvalue"]])
        except ValueError:
            branch_p = None

        dual_sig = (
            fam_p is not None and fam_p < args.family_pcut
            and branch_p is not None and branch_p <= args.branch_pcut
        )
        out_rows.append(row + [("" if fam_p is None else fam_p), dual_sig])

    if missing:
        print(
            f"WARNING: {missing}/{len(rows)} rows had no matching family-wide p-value "
            "(family excluded before the base run, or a results/outdir mismatch)",
            file=sys.stderr,
        )

    new_header = header + ["family_pvalue", "dual_significant"]
    out_all = os.path.join(
        args.outdir, "cafe", "significant_families",
        "combined_changes_per_node.with_family_pvalue.tsv",
    )
    out_sig = os.path.join(
        args.outdir, "cafe", "significant_families",
        "combined_significant_changes_per_node.dual_cutoff.tsv",
    )

    with open(out_all, "w", newline="") as out:
        w = csv.writer(out, delimiter="\t")
        w.writerow(new_header)
        w.writerows(out_rows)

    sig_rows = [r for r in out_rows if r[-1] is True]
    with open(out_sig, "w", newline="") as out:
        w = csv.writer(out, delimiter="\t")
        w.writerow(new_header)
        w.writerows(sig_rows)

    print(f"Wrote {out_all} ({len(out_rows)} rows)")
    print(
        f"Wrote {out_sig} ({len(sig_rows)} rows: "
        f"family p < {args.family_pcut} AND branch p <= {args.branch_pcut})"
    )


if __name__ == "__main__":
    main()
