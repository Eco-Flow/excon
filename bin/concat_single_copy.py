#!/usr/bin/env python3

"""
Build a concatenated (supermatrix) protein alignment from the single-copy
orthogroup alignments produced by OrthoFinder.

OrthoFinder (run with -M msa) writes one aligned FASTA per orthogroup to
MultipleSequenceAlignments/. Sequence headers are gene IDs, so Orthogroups.tsv
is used to map each gene back to the species it came from.

An orthogroup is used only if it is strictly single-copy and complete: every
species must be represented exactly once. Those orthogroups are concatenated in
a stable order and a RAxML-style partition file is emitted alongside, so
IQ-TREE can fit a separate model per orthogroup.

Pure standard library (no BioPython/ete3) to match the rest of bin/.
"""

import argparse
import os
import re
import sys
from collections import OrderedDict


def read_fasta(path):
    """Read a FASTA file into an OrderedDict of {header_first_token: sequence}."""
    seqs = OrderedDict()
    name = None
    chunks = []
    with open(path) as fh:
        for line in fh:
            line = line.rstrip('\n\r')
            if not line:
                continue
            if line.startswith('>'):
                if name is not None:
                    seqs[name] = ''.join(chunks)
                name = line[1:].strip().split()[0]
                chunks = []
            else:
                chunks.append(line.strip())
    if name is not None:
        seqs[name] = ''.join(chunks)
    return seqs


def gene_to_species(orthogroups_tsv):
    """Map every gene ID to its species, using the column headers of Orthogroups.tsv.

    Only usable as a fallback: gene IDs in Orthogroups.tsv are bare (e.g. 'agat-157')
    and are not necessarily unique across species, so a later species can overwrite an
    earlier one. Prefer resolving the species from the alignment header prefix.
    """
    mapping = {}
    with open(orthogroups_tsv) as fh:
        header = fh.readline().rstrip('\n\r').split('\t')
        # First column is the orthogroup ID; the rest are species.
        species = header[1:]
        for line in fh:
            fields = line.rstrip('\n\r').split('\t')
            for sp, cell in zip(species, fields[1:]):
                if not cell.strip():
                    continue
                for gene in cell.split(','):
                    gene = gene.strip()
                    if gene:
                        mapping.setdefault(gene, sp)
    return mapping, species


def species_prefixes(species):
    """Build the header prefixes OrthoFinder writes for each species.

    OrthoFinder names each sequence in MultipleSequenceAlignments/ after the input
    file it came from, with punctuation replaced by underscores — so a species column
    'Genus_species.clean' appears in the alignment as 'Genus_species_clean_<gene>'.
    Returns a list of (prefix, species) sorted longest-first so that a species whose
    name is a prefix of another cannot shadow it.
    """
    prefixes = []
    for sp in species:
        for variant in {sp, sp.replace('.', '_'), re.sub(r'[^A-Za-z0-9]', '_', sp)}:
            prefixes.append((variant + '_', sp))
    prefixes.sort(key=lambda pair: len(pair[0]), reverse=True)
    return prefixes


def resolve_species(header, prefixes, fallback):
    """Resolve a sequence header to a species, by prefix first then by gene ID."""
    for prefix, sp in prefixes:
        if header.startswith(prefix):
            return sp
    return fallback.get(header)


def main():
    parser = argparse.ArgumentParser(
        description='Concatenate single-copy orthogroup alignments into a supermatrix'
    )
    parser.add_argument('-m', '--msa-dir', required=True,
                        help='OrthoFinder MultipleSequenceAlignments/ directory')
    parser.add_argument('-g', '--orthogroups', required=True,
                        help='OrthoFinder Orthogroups.tsv (used to map genes to species)')
    parser.add_argument('-o', '--out-fasta', required=True,
                        help='Output concatenated FASTA alignment')
    parser.add_argument('-p', '--out-partitions', required=True,
                        help='Output RAxML-style partition file')
    parser.add_argument('--min-orthogroups', type=int, default=1,
                        help='Fail if fewer than this many usable orthogroups are found (default: 1)')

    args = parser.parse_args()

    if not os.path.isdir(args.msa_dir):
        sys.exit(
            "ERROR: MSA directory '%s' not found. OrthoFinder only writes "
            "MultipleSequenceAlignments/ when run with -M msa." % args.msa_dir
        )

    mapping, species = gene_to_species(args.orthogroups)
    if not species:
        sys.exit("ERROR: no species columns found in %s" % args.orthogroups)
    species = sorted(species)
    prefixes = species_prefixes(species)

    alignment_files = sorted(
        os.path.join(args.msa_dir, f)
        for f in os.listdir(args.msa_dir)
        if f.endswith(('.fa', '.faa', '.fasta', '.aln'))
    )
    if not alignment_files:
        sys.exit("ERROR: no alignment files found in %s" % args.msa_dir)

    # Accumulate per-species sequence blocks plus the partition ranges.
    blocks = OrderedDict((sp, []) for sp in species)
    partitions = []
    offset = 0
    skipped_multicopy = 0
    skipped_incomplete = 0
    skipped_ragged = 0

    for path in alignment_files:
        og = os.path.splitext(os.path.basename(path))[0]
        seqs = read_fasta(path)
        if not seqs:
            continue

        # Group this orthogroup's sequences by species.
        by_species = {}
        unmapped = False
        for gene, seq in seqs.items():
            sp = resolve_species(gene, prefixes, mapping)
            if sp is None:
                unmapped = True
                break
            by_species.setdefault(sp, []).append(seq)
        if unmapped:
            skipped_incomplete += 1
            continue

        if any(len(v) > 1 for v in by_species.values()):
            skipped_multicopy += 1
            continue
        if len(by_species) != len(species):
            skipped_incomplete += 1
            continue

        lengths = {len(v[0]) for v in by_species.values()}
        if len(lengths) != 1:
            # Not a rectangular alignment — should not happen, but never
            # silently corrupt the supermatrix coordinates.
            skipped_ragged += 1
            continue
        width = lengths.pop()

        for sp in species:
            blocks[sp].append(by_species[sp][0])
        partitions.append((og, offset + 1, offset + width))
        offset += width

    if len(partitions) < args.min_orthogroups:
        sys.exit(
            "ERROR: only %d usable single-copy orthogroup(s) found across %d species "
            "(minimum required: %d). A supermatrix cannot be built."
            % (len(partitions), len(species), args.min_orthogroups)
        )

    with open(args.out_fasta, 'w') as fh:
        for sp in species:
            fh.write('>%s\n' % sp)
            seq = ''.join(blocks[sp])
            for i in range(0, len(seq), 60):
                fh.write(seq[i:i + 60] + '\n')

    with open(args.out_partitions, 'w') as fh:
        for og, start, end in partitions:
            fh.write('AA, %s = %d-%d\n' % (og, start, end))

    print("Supermatrix built successfully")
    print("   Species:                 %d" % len(species))
    print("   Orthogroups used:        %d" % len(partitions))
    print("   Alignment length:        %d" % offset)
    print("   Skipped (multi-copy):    %d" % skipped_multicopy)
    print("   Skipped (incomplete):    %d" % skipped_incomplete)
    if skipped_ragged:
        print("   Skipped (ragged):        %d" % skipped_ragged)
    print("   Alignment written to:    %s" % args.out_fasta)
    print("   Partitions written to:   %s" % args.out_partitions)


if __name__ == '__main__':
    main()
