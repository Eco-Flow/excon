#!/usr/bin/env python3

"""
Write one FASTA per single-copy orthogroup, with each sequence named after the
species it came from.

Orthogroups.tsv says which gene belongs to which species, and the per-species
proteomes hold the sequences, so taking the sequence from the proteome named by
the column keeps the species assignment exact. Gene IDs are only unique within a
species, so identifying species from the sequence headers alone is not reliable.

An orthogroup is written only if it is strictly single-copy: every species in it
represented exactly once. By default it must also be complete (every species
present); --min-occupancy lowers that to a fraction of the species, and the
species an orthogroup lacks are left out of its FASTA (concat_single_copy.py
fills them with gaps).

Pure standard library (no BioPython/ete3) to match the rest of bin/.
"""

import argparse
import math
import os
import sys
from collections import OrderedDict


def read_fasta(path, wanted=None):
    """Read a FASTA file into an OrderedDict of {header_first_token: sequence}.

    When `wanted` is given, only those sequence IDs are kept. Whole proteomes do not
    need to be held in memory at once — only the genes that belong to a single-copy
    orthogroup are ever used.
    """
    seqs = OrderedDict()
    name = None
    chunks = []
    keep = True
    with open(path) as fh:
        for line in fh:
            line = line.rstrip('\n\r')
            if not line:
                continue
            if line.startswith('>'):
                if name is not None and keep:
                    seqs[name] = ''.join(chunks)
                name = line[1:].strip().split()[0]
                keep = wanted is None or name in wanted
                chunks = []
            elif keep:
                chunks.append(line.strip())
    if name is not None and keep:
        seqs[name] = ''.join(chunks)
    return seqs


def find_proteome(proteome_dir, species):
    """Locate the proteome file for a species column of Orthogroups.tsv."""
    candidates = [species, species + '.fasta', species + '.fa', species + '.faa',
                  species + '.fasta.gz']
    for name in candidates:
        path = os.path.join(proteome_dir, name)
        if os.path.isfile(path):
            return path
    # Fall back to any file whose name starts with the species column.
    for name in sorted(os.listdir(proteome_dir)):
        if name.startswith(species):
            return os.path.join(proteome_dir, name)
    return None


def main():
    parser = argparse.ArgumentParser(
        description='Extract single-copy orthogroups as per-species FASTA files'
    )
    parser.add_argument('-g', '--orthogroups', required=True,
                        help='OrthoFinder Orthogroups.tsv')
    parser.add_argument('-p', '--proteome-dir', required=True,
                        help='Directory of per-species proteome FASTA files')
    parser.add_argument('-o', '--out-dir', required=True,
                        help='Output directory for per-orthogroup FASTA files')
    parser.add_argument('--min-occupancy', type=float, default=1.0,
                        help='Smallest fraction of species an orthogroup must contain, each '
                             'exactly once (default: 1.0, every species)')

    args = parser.parse_args()
    if not 0 < args.min_occupancy <= 1:
        sys.exit("ERROR: --min-occupancy must be in (0, 1], got %s" % args.min_occupancy)

    with open(args.orthogroups) as fh:
        header = fh.readline().rstrip('\n\r').split('\t')
        species = header[1:]
        rows = [line.rstrip('\n\r').split('\t') for line in fh if line.strip()]

    if not species:
        sys.exit("ERROR: no species columns found in %s" % args.orthogroups)

    # Round up, so e.g. 0.9 of 25 species needs 23 present, never 22. At least 4,
    # the fewest taxa that carry any information about an unrooted topology.
    min_species = max(min(4, len(species)),
                      math.ceil(args.min_occupancy * len(species) - 1e-9))

    # Work out the single-copy orthogroups first, so only the genes they contain are
    # read from the proteomes. Reading every proteome in full scales with the total
    # gene count across the analysis, which is far more than is needed here.
    single_copy_rows = []
    skipped_not_single = 0
    skipped_low_occupancy = 0
    wanted = {sp: set() for sp in species}

    for fields in rows:
        og = fields[0]
        cells = fields[1:] + [''] * (len(species) - len(fields[1:]))

        genes = {}
        single_copy = True
        for sp, cell in zip(species, cells):
            ids = [g.strip() for g in cell.split(',') if g.strip()]
            if len(ids) > 1:
                single_copy = False
                break
            if ids:
                genes[sp] = ids[0]
        if not single_copy:
            skipped_not_single += 1
            continue
        if len(genes) < min_species:
            skipped_low_occupancy += 1
            continue
        single_copy_rows.append((og, genes))
        for sp, gene in genes.items():
            wanted[sp].add(gene)

    print("Reading %d proteomes for %d single-copy orthogroups"
          % (len(species), len(single_copy_rows)), flush=True)
    proteomes = {}
    for n, sp in enumerate(species, 1):
        path = find_proteome(args.proteome_dir, sp)
        if path is None:
            sys.exit("ERROR: no proteome file found for species '%s' in %s"
                     % (sp, args.proteome_dir))
        proteomes[sp] = read_fasta(path, wanted[sp])
        print("   [%d/%d] %s" % (n, len(species), os.path.basename(path)), flush=True)

    if not os.path.isdir(args.out_dir):
        os.makedirs(args.out_dir)

    written = 0
    missing_gene = 0

    for og, genes in single_copy_rows:

        records = []
        for sp in species:
            if sp not in genes:
                continue
            seq = proteomes[sp].get(genes[sp])
            if seq is None:
                records = None
                break
            records.append((sp, seq))
        if records is None:
            missing_gene += 1
            continue

        with open(os.path.join(args.out_dir, og + '.fa'), 'w') as out:
            for sp, seq in records:
                out.write('>%s\n' % sp)
                for i in range(0, len(seq), 60):
                    out.write(seq[i:i + 60] + '\n')
        written += 1

    if written == 0:
        sys.exit(
            "ERROR: no single-copy orthogroups found across %d species with at least %d "
            "species present. A supermatrix cannot be built. With many species, try "
            "--iqtree_min_occupancy below 1." % (len(species), min_species)
        )

    print("Single-copy orthogroups extracted")
    print("   Species:                     %d" % len(species))
    print("   Minimum species present:     %d (occupancy %s)" % (min_species, args.min_occupancy))
    print("   Orthogroups written:         %d" % written)
    print("   Skipped (not single-copy):   %d" % skipped_not_single)
    print("   Skipped (too few species):   %d" % skipped_low_occupancy)
    if missing_gene:
        print("   Skipped (gene not in proteome): %d" % missing_gene)
    print("   Output directory:            %s" % args.out_dir)


if __name__ == '__main__':
    main()
