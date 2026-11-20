#!/usr/bin/env python3
"""
The fitted parameters of every CAFE5 model the pipeline compared, in one table:
k, lambda, alpha, the error model's epsilon, the root distribution (and the
Poisson root mean), -lnL and the largest lambda CAFE5 allows for the tree.

Reads the k-sweep's Out_cafe_k<k>/ directories (uniform root), the Poisson run
at the best k (Out_cafe_k<k>_poisson/) and its log (CAFE5 writes the fitted
Poisson mean only to the log), and the error model file.

Usage:
  cafe_model_parameters.py --results Out_cafe_k* [--poisson-log cafe_k3_poisson.log]
      [--error-model Base_error_model.txt] [--best-model best_model.txt --best-k best_k.txt]
      -o cafe_model_parameters.tsv
"""

import argparse
import glob
import os
import re


def results_file(d):
    for name in ("Gamma_results.txt", "Base_results.txt"):
        p = os.path.join(d, name)
        if os.path.isfile(p):
            return p
    hits = glob.glob(os.path.join(d, "*_results.txt"))
    hits = [h for h in hits if not h.endswith(("_family_results.txt", "_clade_results.txt"))]
    return hits[0] if hits else None


def parse_results(path):
    """Values from CAFE5's model::write_vital_statistics (and gamma_model's Alpha)."""
    vals = {}
    if not path:
        return vals
    pats = {
        "neg_lnL": r"Final Likelihood \(-lnL\):\s*(\S+)",
        "lambda": r"^Lambda:\s*(\S+)",
        "epsilon": r"^Epsilon:\s*(\S+)",
        "alpha": r"^Alpha:\s*(\S+)",
        "max_possible_lambda": r"^Maximum possible lambda for this topology:\s*(\S+)",
    }
    with open(path) as fh:
        for line in fh:
            for key, pat in pats.items():
                m = re.search(pat, line)
                if m and key not in vals:
                    vals[key] = m.group(1)
    return vals


def error_model_epsilon(path):
    """CAFE5's error model file: 'maxcnt:', 'cntdiff:' lines, then one row per count.
    Epsilon is the last column, shared by all rows (error_model::get_epsilons)."""
    if not path or not os.path.isfile(path) or os.path.getsize(path) == 0:
        return ""
    eps = set()
    with open(path) as fh:
        for line in fh:
            f = line.split()
            if not f or f[0].endswith(":"):
                continue
            eps.add(f[-1])
    return ",".join(sorted(eps))


def poisson_mean(log_path):
    """Last 'Poisson lambda: <x> &  Score: <y>' line of a CAFE5 log."""
    mean = ""
    if log_path and os.path.isfile(log_path):
        with open(log_path, errors="replace") as fh:
            for line in fh:
                m = re.search(r"Poisson lambda:\s*(\S+)", line)
                if m:
                    mean = m.group(1)
    return mean


def main():
    ap = argparse.ArgumentParser(description=__doc__.split("\n\n")[0])
    ap.add_argument("--results", nargs="+", required=True)
    ap.add_argument("--poisson-log")
    ap.add_argument("--error-model")
    ap.add_argument("--best-model")
    ap.add_argument("--best-k")
    ap.add_argument("--fixed", default="",
                    help="comma-separated parameters given rather than estimated "
                         "(lambda, alpha, error_model, poisson_mean)")
    ap.add_argument("-o", "--output", required=True)
    args = ap.parse_args()

    best_root = open(args.best_model).read().strip() if args.best_model and os.path.isfile(args.best_model) else ""
    best_k = open(args.best_k).read().strip() if args.best_k and os.path.isfile(args.best_k) else ""
    eps_file = error_model_epsilon(args.error_model)

    rows = []
    for d in args.results:
        name = os.path.basename(os.path.normpath(d))
        m = re.match(r"Out_cafe_k(\d+)(_poisson)?$", name)
        if not m:
            continue
        k, poisson = m.group(1), bool(m.group(2))
        v = parse_results(results_file(d))
        root = "poisson" if poisson else "uniform"
        rows.append({
            "directory": name,
            "k": k,
            "root_distribution": root,
            "poisson_mean": poisson_mean(args.poisson_log) if poisson else "",
            "lambda": v.get("lambda", "NA"),
            "alpha": v.get("alpha", "") if int(k) > 1 else "",
            "epsilon": v.get("epsilon", eps_file),
            "neg_lnL": v.get("neg_lnL", "NA"),
            "max_possible_lambda": v.get("max_possible_lambda", ""),
            "selected": "BEST" if (k == best_k and root == best_root) else "",
            "fixed_parameters": ",".join(
                f for f in args.fixed.split(",")
                if f and not (f == "alpha" and int(k) == 1)
                and not (f == "poisson_mean" and not poisson)),
        })

    rows.sort(key=lambda r: (int(r["k"]), r["root_distribution"] != "uniform"))
    cols = ["directory", "k", "root_distribution", "poisson_mean", "lambda", "alpha",
            "epsilon", "neg_lnL", "max_possible_lambda", "selected", "fixed_parameters"]
    with open(args.output, "w") as out:
        out.write("\t".join(cols) + "\n")
        for r in rows:
            out.write("\t".join(r[c] for c in cols) + "\n")
    print("Wrote %s (%d fitted models)" % (args.output, len(rows)))


if __name__ == "__main__":
    main()
