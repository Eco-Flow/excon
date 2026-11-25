#!/usr/bin/env python3
"""
KEGG pathway enrichment of the gene families CAFE5 found significantly expanded
or contracted, per node (species or internal node) and direction.

Each family's pathways are the union of its genes' KEGG_Pathway annotations from
eggNOG-mapper (the map##### ids; eggNOG lists each pathway as both ko##### and
map#####). For every node and direction, the significantly changed families
(foreground, from significant_changes_per_node.tsv) are tested against the
background with a one-sided hypergeometric test, P(X >= k), and the p-values are
adjusted with Benjamini-Hochberg over every pathway in that background. The
background is the families CAFE5 modelled (*_change.tab) - for a species, those
that contain it - as for CAFE GO. Only families with at least one pathway count.

Pathway names come from --pathway-names (a TSV of id<TAB>name, as KEGG's
https://rest.kegg.jp/list/pathway returns), or are downloaded from there when
--download-names is given; without either, or if the download fails, the
name column is left empty.

Usage:
  cafe_kegg_enrichment.py --annotations *.emapper.annotations --orthogroups N0.tsv
      --cafe-results Out_cafe --sig-changes significant_changes_per_node.tsv
      [--pathway-names names.tsv | --download-names] -o kegg_enrichment.tsv
"""

import argparse
import glob
import math
import os
import sys
import urllib.request
from collections import defaultdict

KEGG_LIST_URL = "https://rest.kegg.jp/list/pathway"


def species_key(name):
    return os.path.basename(name).split(".")[0]


def read_gene_pathways(paths):
    """species -> gene -> set of map##### ids."""
    out = {}
    for path in paths:
        sp = species_key(path.replace(".emapper.annotations", ""))
        genes = out.setdefault(sp, {})
        col = 12  # eggNOG-mapper 2.1: query ... KEGG_ko(11) KEGG_Pathway(12)
        with open(path) as fh:
            for line in fh:
                if line.startswith("#query"):
                    header = line[1:].rstrip("\r\n").split("\t")
                    if "KEGG_Pathway" in header:
                        col = header.index("KEGG_Pathway")
                    continue
                if line.startswith("#") or not line.strip():
                    continue
                cols = line.rstrip("\r\n").split("\t")
                if len(cols) <= col or cols[col] in ("-", ""):
                    continue
                maps = {p for p in cols[col].split(",") if p.startswith("map")}
                if maps:
                    genes[cols[0]] = maps
    return out


def read_family_pathways(orthogroups, gene_pw):
    """family -> set of pathways, and family -> set of species with genes in it."""
    fam_pw, fam_sp = defaultdict(set), defaultdict(set)
    with open(orthogroups) as fh:
        header = fh.readline().rstrip("\r\n").split("\t")
        sp_cols = [(i, species_key(h)) for i, h in enumerate(header[1:], 1)
                   if h not in ("HOG", "OG", "Gene Tree Parent Clade")]
        for line in fh:
            cols = line.rstrip("\r\n").split("\t")
            fam = cols[0]
            for i, sp in sp_cols:
                if i >= len(cols) or not cols[i].strip():
                    continue
                fam_sp[fam].add(sp)
                for gene in (g.strip() for g in cols[i].split(",")):
                    fam_pw[fam] |= gene_pw.get(sp, {}).get(gene, set())
    return fam_pw, fam_sp


def modelled_families(results_dir):
    fams = set()
    for path in glob.glob(os.path.join(results_dir, "*_change.tab")):
        with open(path) as fh:
            fh.readline()
            fams |= {l.split("\t", 1)[0] for l in fh if l.strip()}
    if not fams:
        sys.exit("ERROR: no *_change.tab in %s" % results_dir)
    return fams


def read_names(path=None, download=False):
    text = ""
    if path and os.path.isfile(path) and os.path.getsize(path):
        text = open(path).read()
    elif download:
        try:
            with urllib.request.urlopen(KEGG_LIST_URL, timeout=60) as r:
                text = r.read().decode()
        except Exception as e:  # offline nodes etc.: keep going without names
            print("WARNING: could not download KEGG pathway names (%s); the name column "
                  "will be empty. Pass --cafe_kegg_pathway_names to give them." % e, file=sys.stderr)
    names = {}
    for line in text.splitlines():
        f = line.split("\t")
        if len(f) >= 2:
            # "map00010<TAB>Glycolysis / Gluconeogenesis" (older lists: "path:map00010")
            pid = f[0].replace("path:", "")
            if pid.startswith("ko"):
                pid = "map" + pid[2:]
            names[pid] = f[1]
    return names


def log_choose(n, k):
    return math.lgamma(n + 1) - math.lgamma(k + 1) - math.lgamma(n - k + 1)


def hypergeom_sf(k, N, K, n):
    """P(X >= k) for X ~ Hypergeometric(population N, K successes, n draws)."""
    hi = min(K, n)
    if k > hi:
        return 0.0
    denom = log_choose(N, n)
    terms = [log_choose(K, x) + log_choose(N - K, n - x) - denom
             for x in range(max(k, n - (N - K)), hi + 1)]
    m = max(terms)
    return min(1.0, math.exp(m) * sum(math.exp(t - m) for t in terms))


def bh(pvalues):
    m = len(pvalues)
    order = sorted(range(m), key=lambda i: pvalues[i])
    adj, prev = [0.0] * m, 1.0
    for rank in range(m, 0, -1):
        i = order[rank - 1]
        prev = min(prev, pvalues[i] * m / rank)
        adj[i] = prev
    return adj


def main():
    ap = argparse.ArgumentParser(description=__doc__.split("\n\n")[0])
    ap.add_argument("--annotations", nargs="+", required=True)
    ap.add_argument("--orthogroups", required=True)
    ap.add_argument("--cafe-results", required=True)
    ap.add_argument("--sig-changes", required=True)
    ap.add_argument("--pathway-names")
    ap.add_argument("--download-names", action="store_true")
    ap.add_argument("--padj-cutoff", type=float, default=0.05)
    ap.add_argument("-o", "--output", required=True)
    ap.add_argument("--significant", default="kegg_enrichment_significant.tsv")
    args = ap.parse_args()

    gene_pw = read_gene_pathways(args.annotations)
    fam_pw, fam_sp = read_family_pathways(args.orthogroups, gene_pw)
    modelled = modelled_families(args.cafe_results)
    names = read_names(args.pathway_names, args.download_names)

    fg = defaultdict(set)  # (node, direction) -> families
    with open(args.sig_changes) as fh:
        header = fh.readline().rstrip("\r\n").split("\t")
        ix = {h: i for i, h in enumerate(header)}
        for line in fh:
            f = line.rstrip("\r\n").split("\t")
            if len(f) == len(header):
                fg[(f[ix["Node"]], f[ix["direction"]])].add(f[ix["HOG"]])

    species = set().union(*fam_sp.values()) if fam_sp else set()
    cols = ["node", "direction", "pathway", "pathway_name", "fg_with_pathway", "fg_annotated",
            "bg_with_pathway", "bg_annotated", "fold_enrichment", "pvalue", "padj_BH", "families"]
    rows = []
    for (node, direction), fams in sorted(fg.items()):
        bg = {f for f in modelled if node not in species or node in fam_sp.get(f, ())}
        bg_ann = {f for f in bg if fam_pw.get(f)}
        fg_ann = {f for f in fams & bg_ann}
        if not fg_ann:
            continue
        N, n = len(bg_ann), len(fg_ann)
        bg_count = defaultdict(int)
        for f in bg_ann:
            for p in fam_pw[f]:
                bg_count[p] += 1
        tests = []
        for p, K in bg_count.items():
            hits = sorted(f for f in fg_ann if p in fam_pw[f])
            tests.append((p, len(hits), K, hypergeom_sf(len(hits), N, K, n), hits))
        adj = bh([t[3] for t in tests])
        for (p, k, K, pv, hits), pa in zip(tests, adj):
            if k == 0:
                continue
            rows.append([node, direction, p, names.get(p, ""), k, n, K, N,
                         round((k / n) / (K / N), 3), pv, pa, ",".join(hits)])

    rows.sort(key=lambda r: (r[0], r[1], r[9]))
    for path, keep in ((args.output, lambda r: True),
                       (args.significant, lambda r: r[10] <= args.padj_cutoff)):
        with open(path, "w") as out:
            out.write("\t".join(cols) + "\n")
            for r in filter(keep, rows):
                out.write("\t".join("%.3g" % v if isinstance(v, float) and c in ("pvalue", "padj_BH") else str(v)
                                    for v, c in zip(r, cols)) + "\n")
    n_sig = sum(1 for r in rows if r[10] <= args.padj_cutoff)
    print("KEGG enrichment: %d node/direction sets tested, %d pathways with BH-adjusted p <= %s"
          % (len(fg), n_sig, args.padj_cutoff))


if __name__ == "__main__":
    main()
