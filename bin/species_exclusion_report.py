#!/usr/bin/env python3

"""
Report which species are associated with orthogroups being excluded from a CAFE
analysis.

CAFE_PREP drops orthogroups that are empty after subsetting, present in a single
species, or that span too wide a range of copy numbers, recording the reason for
each in hog_filtering_report.tsv. This script joins that report back to the
gene-count table and asks, per species, how often it is the one at zero copies
among the excluded orthogroups.

Two very different situations produce exclusions, and the point of this report is
to tell them apart:

  * A poorly annotated genome misses genes across the board. It shows up as a low
    total gene count and as being absent from far more excluded orthogroups than
    its close relatives.
  * Broad taxon sampling excludes orthogroups for real biological reasons. A
    family restricted to one clade looks absent at the root once distant
    outgroups are included, and its copy-number range widens because the distant
    species contribute zeros. Here the absences are spread across many species
    rather than concentrated in one.

Only the first is a data-quality problem worth fixing.

Pure standard library (no BioPython/pandas) to match the rest of bin/.
"""

import argparse
import csv
import os
import sys
from collections import Counter, OrderedDict


def read_counts(path):
    """Read N0.tsv or Orthogroups.tsv into {orthogroup: [count per species]}."""
    with open(path) as fh:
        header = fh.readline().rstrip('\n\r').split('\t')
        # N0.tsv carries three fixed columns before the species; Orthogroups.tsv one.
        fixed = {"HOG", "OG", "Gene Tree Parent Clade", "Orthogroup"}
        start = next((i for i, h in enumerate(header) if h.split('.')[0] not in fixed), 1)
        species = [h.split('.')[0] for h in header[start:]]

        counts = OrderedDict()
        for line in fh:
            fields = line.rstrip('\n\r').split('\t')
            if not fields or not fields[0]:
                continue
            cells = fields[start:] + [''] * (len(species) - len(fields[start:]))
            counts[fields[0]] = [
                len([g for g in cell.split(',') if g.strip()]) for cell in cells
            ]
    return species, counts


def main():
    parser = argparse.ArgumentParser(
        description='Report which species are associated with excluded orthogroups'
    )
    parser.add_argument('-g', '--orthogroups', required=True,
                        help='Gene-count table (N0.tsv or Orthogroups.tsv)')
    parser.add_argument('-r', '--report', required=True,
                        help='hog_filtering_report.tsv written by CAFE_PREP')
    parser.add_argument('-o', '--out', default=None,
                        help='Optional TSV to write the table to as well as stdout')
    parser.add_argument('--low-fraction', type=float, default=0.8,
                        help='Flag species whose gene count is below this fraction '
                             'of the median (default: 0.8)')

    args = parser.parse_args()
    for path in (args.orthogroups, args.report):
        if not os.path.isfile(path):
            sys.exit("ERROR: file not found: %s" % path)

    species, counts = read_counts(args.orthogroups)
    if not species:
        sys.exit("ERROR: no species columns found in %s" % args.orthogroups)

    reasons = {}
    with open(args.report) as fh:
        reader = csv.DictReader(fh, delimiter='\t')
        if 'HOG' not in (reader.fieldnames or []) or \
           'exclusion_reason' not in (reader.fieldnames or []):
            sys.exit("ERROR: %s does not look like a CAFE_PREP filtering report "
                     "(expected HOG and exclusion_reason columns)" % args.report)
        for row in reader:
            reasons[row['HOG']] = row['exclusion_reason']

    known = [h for h in reasons if h in counts]
    if not known:
        sys.exit("ERROR: no orthogroup identifiers are shared between the two files. "
                 "They are probably from different OrthoFinder runs.")

    retained = [h for h in known if reasons[h] == 'retained']
    excluded = [h for h in known if reasons[h] != 'retained']
    reason_types = sorted({reasons[h] for h in excluded})

    print("Orthogroups: %d retained, %d excluded (of %d in the report)"
          % (len(retained), len(excluded), len(reasons)))
    print("Exclusion reasons: %s"
          % ', '.join('%s=%d' % (r, sum(1 for h in excluded if reasons[h] == r))
                      for r in reason_types) or 'none')
    print()

    rows = []
    for i, sp in enumerate(species):
        total = sum(counts[h][i] for h in counts)
        zero_excluded = sum(1 for h in excluded if counts[h][i] == 0)
        zero_retained = sum(1 for h in retained if counts[h][i] == 0)
        per_reason = Counter(reasons[h] for h in excluded if counts[h][i] == 0)
        rows.append((sp, total, zero_excluded, zero_retained, per_reason))

    totals = sorted(r[1] for r in rows)
    median = totals[len(totals) // 2] if totals else 0

    head = "%-38s %9s %9s %9s" % ("species", "genes", "absent_ex", "absent_ret")
    head += ''.join("%22s" % r[:22] for r in reason_types)
    print(head)
    for sp, total, ze, zr, per_reason in sorted(rows, key=lambda r: r[1]):
        line = "%-38s %9d %9d %9d" % (sp, total, ze, zr)
        line += ''.join("%22d" % per_reason.get(r, 0) for r in reason_types)
        if median and total < args.low_fraction * median:
            line += "  <-- low gene count"
        print(line)

    print("\nmedian gene count: %d" % median)
    if excluded:
        spread = [r[2] for r in rows]
        print("absences among excluded orthogroups: min %d, median %d, max %d"
              % (min(spread), sorted(spread)[len(spread) // 2], max(spread)))
        print("\nAbsences concentrated in one or two species suggest an annotation "
              "problem;\nabsences spread evenly suggest the taxon sampling itself is "
              "the cause.")

    if args.out:
        with open(args.out, 'w') as fh:
            fh.write('\t'.join(['species', 'total_genes', 'absent_in_excluded',
                                'absent_in_retained'] + reason_types) + '\n')
            for sp, total, ze, zr, per_reason in sorted(rows, key=lambda r: r[1]):
                fh.write('\t'.join([sp, str(total), str(ze), str(zr)] +
                                   [str(per_reason.get(r, 0)) for r in reason_types]) + '\n')
        print("\nWritten to %s" % args.out)


if __name__ == '__main__':
    main()
