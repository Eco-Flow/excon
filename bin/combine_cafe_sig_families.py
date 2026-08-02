#!/usr/bin/env python3

"""
Concatenate the main-model and large-family significant-families tables into
one combined report spanning every orthogroup CAFE5 was able to fit — the
bulk of families under the shared k/lambda model, plus the high-differential
families that CAFE_RUN_LARGE fits individually (each under its own lambda,
since they don't share one rate with anything else). A "source" column
records which track produced each row so the two remain distinguishable.

Pure standard library, to match the rest of bin/.
"""

import argparse
import csv
import os
import sys


def read_rows(path):
    if not os.path.isfile(path):
        return None, []
    with open(path, newline="") as fh:
        reader = csv.reader(fh, delimiter="\t")
        rows = list(reader)
    if not rows:
        return None, []
    return rows[0], rows[1:]


def combine(main_path, large_path, out_path, source_main, source_large):
    header, main_rows = read_rows(main_path)
    large_header, large_rows = read_rows(large_path)
    header = header or large_header
    if header is None:
        sys.exit("ERROR: neither %s nor %s could be read" % (main_path, large_path))

    with open(out_path, "w", newline="") as out:
        writer = csv.writer(out, delimiter="\t")
        writer.writerow(header + ["source"])
        for row in main_rows:
            writer.writerow(row + [source_main])
        for row in large_rows:
            writer.writerow(row + [source_large])

    return len(main_rows), len(large_rows)


def main():
    ap = argparse.ArgumentParser(description=__doc__)
    ap.add_argument("--main-changes", required=True)
    ap.add_argument("--large-changes", required=True)
    ap.add_argument("--main-sig-changes", required=True)
    ap.add_argument("--large-sig-changes", required=True)
    ap.add_argument("--out-changes", required=True)
    ap.add_argument("--out-sig-changes", required=True)
    args = ap.parse_args()

    n_main, n_large = combine(
        args.main_changes, args.large_changes, args.out_changes,
        "main_model", "large_family_own_lambda",
    )
    print("combined_changes_per_node.tsv: %d main-model rows, %d large-family rows"
          % (n_main, n_large))

    n_main_sig, n_large_sig = combine(
        args.main_sig_changes, args.large_sig_changes, args.out_sig_changes,
        "main_model", "large_family_own_lambda",
    )
    print("combined_significant_changes_per_node.tsv: %d main-model rows, %d large-family rows"
          % (n_main_sig, n_large_sig))


if __name__ == "__main__":
    main()
