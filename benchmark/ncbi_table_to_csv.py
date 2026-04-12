#!/usr/bin/env python3
"""
ncbi_table_to_csv.py
Convert the TSV downloaded from the NCBI Assembly web interface into the
two-column samplesheet format expected by EXCON:
    Species_name,GCF_XXXXXXXXX.X

The NCBI table lists each assembly twice (once as GCA_, once as GCF_).
This script keeps only GCF_ (RefSeq) entries.

By default (--one-per-species), only one assembly is kept per species
(genus + species binomial), choosing the highest GCF number as a proxy
for the most recent assembly. Use --all to keep every strain.

N50 filtering (--min-n50 / --max-n50) lets you target a quality band:
    Contiguous:  --min-n50 500000   (chromosome-level, ≥500 kb)
    Fragmented:  --max-n50 100000   (scaffold assemblies, ≤100 kb)
    Fragmented with floor: --min-n50 10000 --max-n50 100000  (not completely broken)

Usage:
    # One representative per species, contiguous only
    python3 benchmark/ncbi_table_to_csv.py ncbi_table.tsv --min-n50 500000

    # Fragmented band with a quality floor
    python3 benchmark/ncbi_table_to_csv.py ncbi_table.tsv --min-n50 10000 --max-n50 100000

    # All strains in a band
    python3 benchmark/ncbi_table_to_csv.py ncbi_table.tsv --all --max-n50 100000

    # Write directly to a samplesheet
    python3 benchmark/ncbi_table_to_csv.py ncbi_table.tsv \
        --output benchmark/inputs/bacteria_close_contiguous.csv

    # Append to an existing samplesheet
    python3 benchmark/ncbi_table_to_csv.py ncbi_table.tsv \
        --append --output benchmark/inputs/bacteria_close_contiguous.csv
"""

import argparse
import csv
import re
import sys


def organism_to_species_id(name: str) -> str:
    """'Streptococcus pneumoniae TIGR4' -> 'Streptococcus_pneumoniae_TIGR4'"""
    name = name.strip()
    name = re.sub(r'[^A-Za-z0-9]+', '_', name)
    return name.strip('_')


def binomial(name: str) -> str:
    """'Streptococcus pneumoniae TIGR4' -> 'Streptococcus_pneumoniae'"""
    parts = name.strip().split()
    return '_'.join(parts[:2]) if len(parts) >= 2 else name.strip()


def accession_number(acc: str) -> int:
    """GCF_019048645.1 -> 19048645  (for numeric comparison)"""
    m = re.search(r'GCF_0*(\d+)', acc)
    return int(m.group(1)) if m else 0


def main():
    parser = argparse.ArgumentParser(description=__doc__,
                                     formatter_class=argparse.RawDescriptionHelpFormatter)
    parser.add_argument('input', help='NCBI Assembly TSV file')
    parser.add_argument('--output', '-o', default='-',
                        help='Output CSV file (default: stdout)')
    parser.add_argument('--append', '-a', action='store_true',
                        help='Append to output file instead of overwriting')
    parser.add_argument('--all', dest='all_strains', action='store_true',
                        help='Keep all strains; default is one per species (highest GCF number)')
    parser.add_argument('--min-n50', type=int, default=None, metavar='BP',
                        help='Exclude assemblies with Scaffold N50 below this value (e.g. 500000)')
    parser.add_argument('--max-n50', type=int, default=None, metavar='BP',
                        help='Exclude assemblies with Scaffold N50 above this value (e.g. 100000)')
    args = parser.parse_args()

    # --- Read TSV ---
    rows = []
    with open(args.input, newline='', encoding='utf-8-sig') as fh:
        reader = csv.DictReader(fh, delimiter='\t')
        reader.fieldnames = [f.strip() for f in reader.fieldnames]
        for row in reader:
            rows.append({k.strip(): v.strip() for k, v in row.items()})

    if not rows:
        print("ERROR: Input file is empty.", file=sys.stderr)
        sys.exit(1)

    headers = list(rows[0].keys())

    accession_col = next(
        (h for h in headers if 'accession' in h.lower() and 'paired' not in h.lower()), None
    )
    organism_col = next(
        (h for h in headers if 'organism' in h.lower()), None
    )
    n50_col = next(
        (h for h in headers if 'n50' in h.lower()), None
    )

    if not accession_col or not organism_col:
        print(f"ERROR: Could not find required columns in:\n  {headers}", file=sys.stderr)
        print("Expected columns containing 'accession' and 'organism'.", file=sys.stderr)
        sys.exit(1)

    if (args.min_n50 or args.max_n50) and not n50_col:
        print(f"ERROR: --min-n50/--max-n50 requested but no N50 column found in:\n  {headers}",
              file=sys.stderr)
        sys.exit(1)

    if n50_col:
        print(f"N50 column: '{n50_col}'", file=sys.stderr)

    # --- Filter to GCF only, apply N50 band ---
    gcf_rows = []
    seen_acc = set()
    n50_skipped = 0
    for row in rows:
        acc = row.get(accession_col, '').strip()
        org = row.get(organism_col, '').strip()
        if not acc.startswith('GCF_') or acc in seen_acc:
            continue
        seen_acc.add(acc)

        # N50 filter
        if n50_col and (args.min_n50 or args.max_n50):
            raw_n50 = row.get(n50_col, '').strip()
            try:
                n50 = int(raw_n50)
            except ValueError:
                n50_skipped += 1
                continue  # skip rows with missing N50
            if args.min_n50 and n50 < args.min_n50:
                n50_skipped += 1
                continue
            if args.max_n50 and n50 > args.max_n50:
                n50_skipped += 1
                continue

        gcf_rows.append((org, acc))

    if n50_skipped:
        print(f"Skipped {n50_skipped} assemblies outside N50 range "
              f"(min={args.min_n50}, max={args.max_n50}).", file=sys.stderr)

    if not gcf_rows:
        print("WARNING: No GCF_ entries found.", file=sys.stderr)
        sys.exit(1)

    # --- One per species: keep highest GCF number per binomial ---
    if not args.all_strains:
        best: dict[str, tuple[str, str, int]] = {}  # binomial -> (org, acc, number)
        for org, acc in gcf_rows:
            key = binomial(org)
            num = accession_number(acc)
            if key not in best or num > best[key][2]:
                best[key] = (org, acc, num)
        gcf_rows = [(org, acc) for org, acc, _ in best.values()]
        print(f"Kept {len(gcf_rows)} species (one representative each, highest GCF number).",
              file=sys.stderr)
    else:
        print(f"Kept {len(gcf_rows)} assemblies (all strains).", file=sys.stderr)

    entries = [(organism_to_species_id(org), acc) for org, acc in gcf_rows]

    # --- Write output ---
    if args.output == '-':
        for species_id, acc in entries:
            print(f"{species_id},{acc}")
    else:
        mode = 'a' if args.append else 'w'
        with open(args.output, mode, newline='') as fh:
            for species_id, acc in entries:
                fh.write(f"{species_id},{acc}\n")
        action = "Appended" if args.append else "Written"
        print(f"{action} {len(entries)} entries to {args.output}", file=sys.stderr)


if __name__ == '__main__':
    main()
