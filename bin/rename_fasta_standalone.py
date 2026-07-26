#!/usr/bin/env python3

"""
Rebuild the *.clean.fasta proteomes outside of Nextflow, from the published
GFFREAD proteins and AGAT annotations.

This is the same transcript-to-gene renaming that the RENAME_FASTA process
applies, split out so the proteomes can be recovered when a work directory has
been deleted but results/gffread/ and results/agat/ survive. The gene IDs it
produces are the ones that appear in Orthogroups.tsv/N0.tsv, so the output is
what --proteome_dir expects.

Pure standard library (no BioPython/ete3) to match the rest of bin/.
"""

import argparse
import os
import re
import sys


def transcript_to_gene(gff_path):
    """Map transcript IDs to gene IDs, covering NCBI, Ensembl and Braker styles."""
    mapping = {}
    with open(gff_path) as fh:
        for line in fh:
            if line.startswith('#'):
                continue
            parts = line.rstrip('\n\r').split('\t')
            if len(parts) < 9 or parts[2] != 'mRNA':
                continue
            attrs = {}
            for a in parts[8].split(';'):
                if '=' in a:
                    k, v = a.split('=', 1)
                    attrs[k.strip()] = v.strip()

            tran_id = attrs.get('ID', '')
            gene_id = attrs.get('Parent', '')

            if ':' in gene_id:
                gene_id = gene_id.split(':')[-1]
            gene_id = gene_id.replace('gene-', '')

            tran_id = tran_id.replace('transcript:', '')
            gene_id = gene_id.replace('gene:', '')

            if not gene_id and '.' in tran_id:
                gene_id = tran_id.rsplit('.', 1)[0]

            if tran_id and gene_id:
                mapping[tran_id] = gene_id
                mapping['rna-' + tran_id] = gene_id
                mapping['transcript:' + tran_id] = gene_id
                mapping[tran_id.replace('rna-', '')] = gene_id
    return mapping


def rename(fasta_path, mapping, out_path):
    seen = set()
    written = 0
    duplicates = 0
    with open(fasta_path) as fin, open(out_path, 'w') as fout:
        skip = False
        for line in fin:
            if line.startswith('>'):
                seq_id = line[1:].strip().split()[0]
                gene_id = mapping.get(seq_id, seq_id)
                if gene_id == seq_id and '.' in seq_id:
                    gene_id = seq_id.rsplit('.', 1)[0]
                if gene_id in seen:
                    skip = True
                    duplicates += 1
                else:
                    seen.add(gene_id)
                    skip = False
                    fout.write('>%s\n' % gene_id)
                    written += 1
            else:
                if not skip:
                    fout.write(line.replace('.', '').replace('*', ''))
    return written, duplicates


def strip_ext(name):
    return re.sub(r'\.(fa|faa|fasta|gff|gff3)(\.gz)?$', '', name)


def main():
    parser = argparse.ArgumentParser(
        description='Rebuild *.clean.fasta proteomes from published GFFREAD and AGAT output'
    )
    parser.add_argument('-f', '--gffread-dir', required=True,
                        help='Directory of GFFREAD protein FASTA files (results/gffread)')
    parser.add_argument('-g', '--agat-dir', required=True,
                        help='Directory of AGAT GFF files (results/agat)')
    parser.add_argument('-o', '--out-dir', required=True,
                        help='Directory to write <species>.clean.fasta into')

    args = parser.parse_args()

    fastas = {strip_ext(f): os.path.join(args.gffread_dir, f)
              for f in sorted(os.listdir(args.gffread_dir))
              if f.endswith(('.fa', '.faa', '.fasta'))}
    gffs = {strip_ext(f): os.path.join(args.agat_dir, f)
            for f in sorted(os.listdir(args.agat_dir))
            if f.endswith(('.gff', '.gff3'))}

    if not fastas:
        sys.exit("ERROR: no FASTA files found in %s" % args.gffread_dir)

    if not os.path.isdir(args.out_dir):
        os.makedirs(args.out_dir)

    missing = []
    total = 0
    for species, fasta in fastas.items():
        gff = gffs.get(species)
        if gff is None:
            # AGAT filenames often carry a suffix; fall back to a prefix match.
            candidates = [v for k, v in gffs.items() if k.startswith(species) or species.startswith(k)]
            gff = candidates[0] if len(candidates) == 1 else None
        if gff is None:
            missing.append(species)
            continue

        mapping = transcript_to_gene(gff)
        out = os.path.join(args.out_dir, species + '.clean.fasta')
        written, duplicates = rename(fasta, mapping, out)
        total += written
        note = '  (%d duplicate gene IDs skipped)' % duplicates if duplicates else ''
        print("   %-45s %6d sequences%s" % (species + '.clean.fasta', written, note))

    if missing:
        sys.exit("\nERROR: no matching GFF found for: %s" % ', '.join(missing))

    print("Rebuilt %d proteomes (%d sequences) into %s" % (len(fastas), total, args.out_dir))


if __name__ == '__main__':
    main()
