#!/usr/bin/env python3
"""
Find gene families made mostly of transposable elements (TEs), from eggNOG-mapper
annotations, and remove them from the orthogroup table before CAFE5.

TE families can dominate expansions and contractions, and are rarely the "genes"
of interest. A gene counts as a TE when any of its Pfam domains, GO terms or its
eggNOG description matches a signature in the signatures file
(assets/te_signatures.tsv by default). A family is removed when at least
--min-fraction of its eggNOG-annotated genes are TEs; families with no annotated
genes are kept, since nothing is known about them.

Writes:
  <out-dir>/<orthogroups file name>  the orthogroup table without the TE families
  te_families.tsv                    every family with at least one TE gene: its
                                     counts, fraction, whether it was removed, and
                                     the signatures that matched

Usage:
  te_family_filter.py --orthogroups Orthogroups.tsv --annotations *.emapper.annotations
      [--signatures te_signatures.tsv] [--min-fraction 0.5] --out-dir filtered
"""

import argparse
import os
import re
import sys
from collections import Counter


def read_signatures(path):
    sig = {"pfam": [], "go": set(), "description": []}
    with open(path) as fh:
        for line in fh:
            if not line.strip() or line.startswith("#"):
                continue
            f = line.rstrip("\r\n").split("\t")
            kind, pattern = f[0].strip(), f[1].strip() if len(f) > 1 else ""
            if not pattern:
                continue
            if kind == "pfam":
                sig["pfam"].append(re.compile(r"^(?:%s)$" % pattern))
            elif kind == "go":
                sig["go"].add(pattern)
            elif kind == "description":
                sig["description"].append(re.compile(pattern, re.I))
            else:
                sys.exit("ERROR: unknown signature type '%s' in %s (use pfam, go or description)"
                         % (kind, path))
    return sig


def species_key(name):
    """Species name as written in both file names and table columns: up to the first '.'"""
    return os.path.basename(name).split(".")[0]


def read_annotations(paths):
    """species -> gene -> (pfams, gos, description) from *.emapper.annotations."""
    annot = {}
    for path in paths:
        sp = species_key(path.replace(".emapper.annotations", ""))
        genes = annot.setdefault(sp, {})
        header = None
        with open(path) as fh:
            for line in fh:
                if line.startswith("#query"):
                    header = line[1:].rstrip("\r\n").split("\t")
                    continue
                if line.startswith("#") or not line.strip():
                    continue
                cols = line.rstrip("\r\n").split("\t")
                if header:
                    row = dict(zip(header, cols))
                    desc, gos, pfams = row.get("Description", "-"), row.get("GOs", "-"), row.get("PFAMs", "-")
                else:  # emapper 2.1 column order
                    desc = cols[7] if len(cols) > 7 else "-"
                    gos = cols[9] if len(cols) > 9 else "-"
                    pfams = cols[-1]
                genes[cols[0]] = (
                    [] if pfams in ("-", "") else pfams.split(","),
                    set() if gos in ("-", "") else set(gos.split(",")),
                    "" if desc == "-" else desc,
                )
    return annot


def te_evidence(info, sig):
    pfams, gos, desc = info
    hits = ["pfam:" + p for p in pfams if any(r.match(p) for r in sig["pfam"])]
    hits += ["go:" + g for g in sorted(gos & sig["go"])]
    if desc and any(r.search(desc) for r in sig["description"]):
        hits.append("description")
    return hits


def main():
    ap = argparse.ArgumentParser(description=__doc__.split("\n\n")[0])
    ap.add_argument("--orthogroups", required=True)
    ap.add_argument("--annotations", nargs="+", required=True)
    ap.add_argument("--signatures", required=True)
    ap.add_argument("--min-fraction", type=float, default=0.5)
    ap.add_argument("--out-dir", required=True)
    ap.add_argument("--report", default="te_families.tsv")
    args = ap.parse_args()
    if not 0 < args.min_fraction <= 1:
        sys.exit("ERROR: --min-fraction must be in (0, 1], got %s" % args.min_fraction)

    sig = read_signatures(args.signatures)
    annot = read_annotations(args.annotations)
    if not annot:
        sys.exit("ERROR: no eggNOG annotations read")

    with open(args.orthogroups) as fh:
        lines = fh.read().splitlines()
    header = lines[0].split("\t")
    # Species columns: those named like an annotation file. N0.tsv (OrthoFinder 2) has
    # HOG, OG and Gene Tree Parent Clade before them; Orthogroups.tsv only the OG id.
    sp_cols = [(i, species_key(h)) for i, h in enumerate(header) if species_key(h) in annot]
    if not sp_cols:
        sys.exit("ERROR: no column of %s matches an annotated species (%s)"
                 % (args.orthogroups, ", ".join(sorted(annot)[:5])))

    os.makedirs(args.out_dir, exist_ok=True)
    kept, removed = [lines[0]], 0
    report = []
    for line in lines[1:]:
        if not line.strip():
            continue
        cols = line.split("\t")
        n_genes = n_annot = n_te = 0
        evidence = Counter()
        for i, sp in sp_cols:
            if i >= len(cols) or not cols[i].strip():
                continue
            for gene in (g.strip() for g in cols[i].split(",")):
                if not gene:
                    continue
                n_genes += 1
                info = annot[sp].get(gene)
                if info is None:
                    continue
                n_annot += 1
                hits = te_evidence(info, sig)
                if hits:
                    n_te += 1
                    evidence.update(hits)
        frac = n_te / n_annot if n_annot else 0.0
        is_te = n_te > 0 and frac >= args.min_fraction
        if n_te:
            report.append((cols[0], n_genes, n_annot, n_te, frac, is_te,
                           ",".join("%s(%d)" % (k, v) for k, v in evidence.most_common(5))))
        if is_te:
            removed += 1
        else:
            kept.append(line)

    out_table = os.path.join(args.out_dir, os.path.basename(args.orthogroups))
    with open(out_table, "w") as out:
        out.write("\n".join(kept) + "\n")
    with open(args.report, "w") as out:
        out.write("HOG\tn_genes\tn_annotated\tn_te_genes\tte_fraction\tremoved\tevidence\n")
        for hog, n, na, nt, fr, te, ev in sorted(report, key=lambda r: (-r[4], r[0])):
            out.write("%s\t%d\t%d\t%d\t%.3f\t%s\t%s\n" % (hog, n, na, nt, fr, "yes" if te else "no", ev))

    total = len(lines) - 1
    print("TE filter: removed %d of %d families before CAFE5 (at least %.0f%% of their "
          "eggNOG-annotated genes are transposable elements); %d more have some TE genes and "
          "were kept. See cafe/te_filter/te_families.tsv."
          % (removed, total, 100 * args.min_fraction, len(report) - removed, ))


if __name__ == "__main__":
    main()
