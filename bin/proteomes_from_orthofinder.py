#!/usr/bin/env python3

"""
Reconstruct the per-species proteomes that were given to OrthoFinder, from a
completed OrthoFinder run.

OrthoFinder keeps the input sequences in WorkingDirectory/SpeciesN.fa under
internal identifiers, together with SpeciesIDs.txt (index -> input filename) and
SequenceIDs.txt (internal id -> original sequence header). Together these restore
the original proteomes exactly, which is useful when the pipeline work directory
has been deleted but the published OrthoFinder results remain.

The output is what --proteome_dir expects.

Pure standard library (no BioPython/ete3) to match the rest of bin/.
"""

import argparse
import os
import sys


def read_map(path, sep=': '):
    """Read an OrthoFinder 'key: value' index file."""
    mapping = {}
    with open(path) as fh:
        for line in fh:
            line = line.rstrip('\n\r')
            if not line:
                continue
            key, _, value = line.partition(sep)
            if not _:
                continue
            mapping[key.strip()] = value.strip()
    return mapping


def main():
    parser = argparse.ArgumentParser(
        description='Recover per-species proteomes from a completed OrthoFinder run'
    )
    parser.add_argument('-r', '--orthofinder-results', required=True,
                        help='OrthoFinder results directory (containing WorkingDirectory/)')
    parser.add_argument('-o', '--out-dir', required=True,
                        help='Directory to write the per-species FASTA files into')
    parser.add_argument('--line-width', type=int, default=0,
                        help='Wrap sequences at this width (0 = one line per sequence)')

    args = parser.parse_args()

    wd = args.orthofinder_results
    if not os.path.isfile(os.path.join(wd, 'SpeciesIDs.txt')):
        wd = os.path.join(args.orthofinder_results, 'WorkingDirectory')
    for required in ('SpeciesIDs.txt', 'SequenceIDs.txt'):
        if not os.path.isfile(os.path.join(wd, required)):
            sys.exit("ERROR: %s not found under '%s'. Point --orthofinder-results at an "
                     "OrthoFinder results directory containing WorkingDirectory/."
                     % (required, args.orthofinder_results))

    species = read_map(os.path.join(wd, 'SpeciesIDs.txt'))
    sequences = read_map(os.path.join(wd, 'SequenceIDs.txt'))

    if not os.path.isdir(args.out_dir):
        os.makedirs(args.out_dir)

    total = 0
    for index, filename in sorted(species.items(), key=lambda kv: int(kv[0])):
        src = os.path.join(wd, 'Species%s.fa' % index)
        if not os.path.isfile(src):
            sys.exit("ERROR: %s not found; the OrthoFinder run looks incomplete." % src)

        written = 0
        with open(src) as fin, open(os.path.join(args.out_dir, filename), 'w') as fout:
            for line in fin:
                if line.startswith('>'):
                    internal = line[1:].strip().split()[0]
                    original = sequences.get(internal)
                    if original is None:
                        sys.exit("ERROR: sequence id '%s' is missing from SequenceIDs.txt"
                                 % internal)
                    # OrthoFinder stores the full original header; keep the first token,
                    # which is the gene ID used throughout Orthogroups.tsv.
                    fout.write('>%s\n' % original.split()[0])
                    written += 1
                else:
                    fout.write(line)
        total += written
        print("   %-45s %6d sequences" % (filename, written))

    print("Recovered %d proteomes (%d sequences) into %s"
          % (len(species), total, args.out_dir))


if __name__ == '__main__':
    main()
