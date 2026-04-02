#!/usr/bin/env python3
import os
import re
import sys
import math
import platform

def parse_score(result_dir):
    for fname in ["Base_results.txt", "Gamma_results.txt"]:
        fpath = os.path.join(result_dir, fname)
        if os.path.isfile(fpath):
            with open(fpath) as fh:
                for line in fh:
                    m = re.match(r"Model \w+ Final Likelihood \(-lnL\):\s+([\d.]+)", line)
                    if m:
                        return float(m.group(1))
    return None

def parse_n_families(result_dir):
    for fname in ["Base_family_results.txt", "Gamma_family_results.txt"]:
        fpath = os.path.join(result_dir, fname)
        if os.path.isfile(fpath):
            with open(fpath) as fh:
                lines = [l for l in fh if not l.startswith("#") and l.strip()]
            return len(lines)
    return None

dirs = sorted(
    [d for d in os.listdir(".") if re.match(r"Out_cafe_k\d+", d) and os.path.isdir(d)],
    key=lambda d: int(re.search(r"(\d+)", d).group(1))
)

if not dirs:
    sys.exit("ERROR: No Out_cafe_k* directories found in working directory")

rows = []
for d in dirs:
    k     = int(re.search(r"(\d+)", d).group(1))
    score = parse_score(d)

    if score is None:
        print(f"WARNING: no score found in {d} — skipping", file=sys.stderr)
        continue

    n_params   = 1 if k == 1 else 2
    aic        = 2 * n_params + 2 * score
    n_families = parse_n_families(d)
    bic        = (n_params * math.log(n_families) + 2 * score) if n_families else None

    rows.append({
        "k":          k,
        "n_params":   n_params,
        "neg_lnL":    score,
        "AIC":        aic,
        "BIC":        bic,
        "n_families": n_families,
    })

if not rows:
    sys.exit("ERROR: Could not extract scores from any CAFE5 results directory")

best_row = min(rows, key=lambda r: r["AIC"])
best_k   = best_row["k"]

with open("model_selection.tsv", "w") as fh:
    fh.write("k\tn_params\tneg_lnL\tAIC\tBIC\tn_families\tSelected\n")
    for r in sorted(rows, key=lambda r: r["k"]):
        selected = "BEST" if r["k"] == best_k else ""
        bic_str  = f"{r['BIC']:.4f}" if r["BIC"] is not None else "NA"
        fam_str  = str(r["n_families"]) if r["n_families"] is not None else "NA"
        fh.write(
            f"{r['k']}\t{r['n_params']}\t{r['neg_lnL']:.4f}\t"
            f"{r['AIC']:.4f}\t{bic_str}\t{fam_str}\t{selected}\n"
        )

with open("best_k.txt", "w") as fh:
    fh.write(str(best_k) + "\n")

print(f"Model selection complete — best k = {best_k} (AIC = {best_row['AIC']:.4f})")

with open("versions.yml", "w") as fh:
    fh.write("CAFE_SELECT_K:\n")
    fh.write(f"    python: {platform.python_version()}\n")
