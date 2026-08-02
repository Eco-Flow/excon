#!/usr/bin/env python3

"""
Split hog_gene_counts_large.tsv (one row per high-differential family) into
one TSV per family, so CAFE_RUN_LARGE can run each independently and fit its
own lambda instead of sharing one across every excluded family at once.

Pure standard library, to match the rest of bin/.
"""

import argparse
import os
import re
import sys

SAFE = re.compile(r"[^A-Za-z0-9_.-]")


def main():
    ap = argparse.ArgumentParser(description=__doc__)
    ap.add_argument("large_counts", help="hog_gene_counts_large.tsv")
    ap.add_argument("-o", "--outdir", default="large_family_splits")
    args = ap.parse_args()

    os.makedirs(args.outdir, exist_ok=True)

    with open(args.large_counts) as fh:
        header = fh.readline()
        if not header:
            sys.exit("ERROR: %s is empty" % args.large_counts)

        n = 0
        for line in fh:
            if not line.strip():
                continue
            hog = line.rstrip("\n").split("\t")[1]
            safe_hog = SAFE.sub("_", hog)
            with open(os.path.join(args.outdir, safe_hog + ".tsv"), "w") as out:
                out.write(header)
                out.write(line)
            n += 1

    print("Split %d families into %s/" % (n, args.outdir))


if __name__ == "__main__":
    main()
