#!/usr/bin/env python3
"""
One row per gene family saying what became of it in the CAFE5 analysis.

Until now this was spread over hog_filtering_report.tsv (left out before CAFE5),
the CAFE5 log (dropped at the root), CAFE_RUN_LARGE's stderr (did not converge)
and nothing at all (large-family tasks that failed and were ignored).

Usage:
  cafe_family_status.py --filtering-report hog_filtering_report.tsv \
      --best-results best_cafe_results \
      [--large-counts hog_gene_counts_large.tsv] [--large-status large_status_*.tsv ...] \
      -o cafe_family_status.tsv

status is one of:
  modelled                     in the main CAFE5 model (the selected best model)
  modelled_large               fitted by the large-family track
  excluded_single_species      genes in one species only; nothing to model
  excluded_empty_after_subset  no genes in the species analysed (--cafe_clade/--cafe_species)
  excluded_transposable_element  mostly transposable elements (--cafe_te_filter), removed
                               before CAFE_PREP
  not_present_at_root          no genes on one side of the root, so CAFE5 dropped it
                               from the main model (--cafe_zero_root fits these in a
                               separate sensitivity run, cafe/zero_root/)
  large_not_present_at_root    the same, in the large-family track
  large_did_not_converge       the large-family fit did not reach a finite likelihood
  large_task_failed            the large-family task failed on every attempt
  not_in_results               given to CAFE5 but missing from its results for another reason
"""

import argparse
import glob
import os
import sys
from collections import Counter, OrderedDict

def is_large_reason(reason):
    """Left out of the main run for its size: max_copies_ge_<--cafe_max_copies> or
    differential_gt_threshold."""
    return reason.startswith("max_copies_ge_") or reason == "differential_gt_threshold"


def read_tsv(path):
    with open(path) as fh:
        header = fh.readline().rstrip("\r\n").split("\t")
        for line in fh:
            if line.strip():
                yield dict(zip(header, line.rstrip("\r\n").split("\t")))


def change_tab_families(results_dir):
    """Families in a CAFE5 results directory: the first column of *_change.tab."""
    found = set()
    paths = glob.glob(os.path.join(results_dir, "*_change.tab"))
    if not paths:
        sys.exit("ERROR: no *_change.tab in %s" % results_dir)
    for path in paths:
        with open(path) as fh:
            fh.readline()
            for line in fh:
                if line.strip():
                    found.add(line.split("\t", 1)[0])
    return found


def main():
    ap = argparse.ArgumentParser(description=__doc__.split("\n\n")[0])
    ap.add_argument("--filtering-report", required=True)
    ap.add_argument("--best-results", required=True,
                    help="CAFE5 results directory of the selected main model")
    ap.add_argument("--large-counts",
                    help="hog_gene_counts_large.tsv: the families sent to the large-family track")
    ap.add_argument("--large-status", nargs="*", default=[],
                    help="large_status_*.tsv files written by CAFE_RUN_LARGE")
    ap.add_argument("--te-report",
                    help="te_families.tsv from TE_FAMILY_FILTER (--cafe_te_filter)")
    ap.add_argument("-o", "--output", required=True)
    ap.add_argument("--summary", default="cafe_family_status_summary.txt")
    args = ap.parse_args()

    report = OrderedDict((row["HOG"], row) for row in read_tsv(args.filtering_report))
    modelled = change_tab_families(args.best_results)

    large_sent = set()
    if args.large_counts and os.path.isfile(args.large_counts) and os.path.getsize(args.large_counts):
        large_sent = {row["HOG"] for row in read_tsv(args.large_counts)}

    large_outcome = {}
    for path in args.large_status:
        with open(path) as fh:
            for line in fh:
                if line.strip():
                    hog, outcome = line.rstrip("\r\n").split("\t")[:2]
                    large_outcome[hog] = outcome

    rows = []
    for hog, rec in report.items():
        reason = rec.get("exclusion_reason", "")
        if reason == "retained":
            track = "main"
            if hog in modelled:
                status = "modelled"
            else:
                # The main model never runs with -z, so a retained family missing
                # from its results had no genes on one side of the root
                status = "not_present_at_root"
        elif is_large_reason(reason):
            track = "large"
            outcome = large_outcome.get(hog)
            if outcome == "fitted":
                status = "modelled_large"
            elif outcome in ("not_present_at_root", "did_not_converge"):
                status = "large_" + outcome
            elif hog in large_sent:
                status = "large_task_failed"
            else:
                # Not sent to the large-family track at all (should not happen)
                status = "not_in_results"
        elif reason in ("single_species", "empty_after_subset"):
            track = "none"
            status = "excluded_" + reason
        else:
            track = "none"
            status = "excluded_" + (reason or "unknown")
        rows.append((hog, status, track, reason, rec))

    # Families removed as transposable elements never reached CAFE_PREP's report
    if args.te_report:
        for row in read_tsv(args.te_report):
            if row.get("removed") == "yes" and row["HOG"] not in report:
                rows.append((row["HOG"], "excluded_transposable_element", "none",
                             "transposable_element", {}))

    cols = ["n_species_present", "n_max", "n_min", "total_genes", "differential"]
    with open(args.output, "w") as out:
        out.write("HOG\tstatus\ttrack\tprep_filter\t" + "\t".join(cols) + "\n")
        for hog, status, track, reason, rec in rows:
            out.write("\t".join([hog, status, track, reason] + [rec.get(c, "") for c in cols]) + "\n")

    counts = Counter(r[1] for r in rows)
    order = ["modelled", "modelled_large", "excluded_single_species", "excluded_empty_after_subset",
             "excluded_transposable_element",
             "not_present_at_root", "large_not_present_at_root", "large_did_not_converge",
             "large_task_failed", "not_in_results"]
    order += sorted(set(counts) - set(order))
    lines = ["CAFE family status: %d families in the orthogroup table" % len(rows)]
    lines += ["  %-28s %d" % (s, counts[s]) for s in order if counts[s]]
    lines.append("Per-family detail: cafe/cafe_family_status.tsv")
    with open(args.summary, "w") as fh:
        fh.write("\n".join(lines) + "\n")
    print("\n".join(lines))


if __name__ == "__main__":
    main()
