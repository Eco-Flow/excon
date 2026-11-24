#!/usr/bin/env python3
"""
Prepare a user-supplied proteome (--input_proteins) for OrthoFinder, in place of
the genome + annotation path (AGAT -> GFFREAD -> RENAME_FASTA).

The proteome should hold one protein per gene (e.g. the longest isoform), named
by the ID to use for that gene. This step:
  * reads plain or gzipped FASTA;
  * treats '.' as a stop codon as well as '*' ('.' is not an amino acid, and
    gffread -y writes stops as '.' unless run with -S);
  * removes a terminal stop, and handles internal stops with
    --internal-stop-action, as RENAME_FASTA does: 'longest_orf' (keep the
    longest stretch between stops), 'strip' (remove the stops) or 'drop' (drop
    the protein);
  * drops empty sequences and second copies of an ID, with a warning;
  * writes <species>.clean.fasta, <species>.internal_stop_codons.tsv when any
    protein had an internal stop, and <species>.protein_ids.gff, a minimal GFF
    mapping each protein to itself, so EGGNOG_TO_GO keeps the IDs as they are.

Usage: prepare_proteins.py <species> <proteins.fa[.gz]> <internal_stop_action>
"""

import gzip
import sys


def read_fasta(path):
    opener = gzip.open if path.endswith(".gz") else open
    name, chunks = None, []
    with opener(path, "rt") as fh:
        for line in fh:
            line = line.strip()
            if not line:
                continue
            if line.startswith(">"):
                if name is not None:
                    yield name, "".join(chunks)
                name, chunks = line[1:].split()[0] if len(line) > 1 else "", []
            else:
                chunks.append(line)
    if name is not None:
        yield name, "".join(chunks)


def main():
    if len(sys.argv) != 4:
        sys.exit(__doc__)
    species, path, action = sys.argv[1:]
    if action not in ("strip", "drop", "longest_orf"):
        sys.exit("ERROR: internal stop action must be strip, drop or longest_orf, got %s" % action)

    seen, stops, n_in, n_out, dupes, empty = set(), [], 0, 0, 0, 0
    with open(species + ".clean.fasta", "w") as fa, open(species + ".protein_ids.gff", "w") as gff:
        gff.write("##gff-version 3\n")
        for name, seq in read_fasta(path):
            n_in += 1
            if not name:
                empty += 1
                continue
            if name in seen:
                dupes += 1
                print("WARNING: %s: second protein named %s skipped" % (species, name), flush=True)
                continue
            seq = seq.upper().replace(".", "*")
            body = seq.rstrip("*")
            if "*" in body:
                original = len(body)
                if action == "drop":
                    stops.append((name, original, 0))
                    continue
                body = max(body.split("*"), key=len) if action == "longest_orf" else body.replace("*", "")
                stops.append((name, original, len(body)))
            if not body:
                empty += 1
                continue
            seen.add(name)
            fa.write(">%s\n%s\n" % (name, body))
            gff.write("%s\tinput_proteins\tmRNA\t1\t%d\t.\t+\t.\tID=%s;Parent=%s\n"
                      % (name, 3 * len(body), name, name))
            n_out += 1

    if stops:
        with open(species + ".internal_stop_codons.tsv", "w") as rep:
            rep.write("gene_id\taction\toriginal_length\tkept_length\n")
            for g, o, k in stops:
                rep.write("%s\t%s\t%d\t%d\n" % (g, action, o, k))
    print("%s: %d proteins read, %d written; %d had an internal stop (%s), %d duplicate IDs "
          "and %d empty sequences skipped" % (species, n_in, n_out, len(stops), action, dupes, empty))
    if n_out == 0:
        sys.exit("ERROR: no proteins left for %s" % species)


if __name__ == "__main__":
    main()
