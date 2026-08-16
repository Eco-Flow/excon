process RENAME_FASTA {
    tag "$meta.id"
    label 'process_single'

    container "${ workflow.containerEngine == 'singularity' && !task.ext.singularity_pull_docker_container ?
        'https://depot.galaxyproject.org/singularity/eggnog-mapper:2.1.13--pyhdfd78af_2' :
        'biocontainers/eggnog-mapper:2.1.13--pyhdfd78af_2' }"

    input:
    tuple val(meta), path(fasta)
    tuple val(meta2), path(gff)

    output:
    tuple val(meta), path("${meta.id}.clean.fasta"), emit: fasta
    tuple val(meta), path("${meta.id}.internal_stop_codons.tsv"), emit: internal_stop_report, optional: true
    tuple val("${task.process}"), val('python'), eval("python3 --version | sed 's/Python //'"), emit: versions, topic: versions

    script:
    """
    python3 <<EOF
    # 'strip' (default) removes every '*' including internal ones, splicing the two
    # flanking peptide fragments together; 'drop' discards the whole gene instead.
    # See --internal_stop_action in nextflow.config for the reasoning.
    INTERNAL_STOP_ACTION = "${params.internal_stop_action}"

    tran_to_gene = {}
    with open("${gff}") as f:
        for line in f:
            if line.startswith("#"): continue
            parts = line.strip().split("\\t")
            if len(parts) < 9 or parts[2] != "mRNA": continue
            attrs = {}
            for a in parts[8].split(";"):
                if "=" in a:
                    k, v = a.split("=", 1)
                    attrs[k.strip()] = v.strip()

            tran_id = attrs.get("ID", "")
            gene_id = attrs.get("Parent", "")

            # Handle NCBI style: Parent=gene-LOC123
            if ":" in gene_id:
                gene_id = gene_id.split(":")[-1]
            gene_id = gene_id.replace("gene-", "")

            # Handle Ensembl style: ID=transcript:ENST, Parent=gene:ENSG
            tran_id = tran_id.replace("transcript:", "")
            gene_id = gene_id.replace("gene:", "")

            # Handle Braker/augustus style: ID=g1.t1 -> gene=g1
            if not gene_id and "." in tran_id:
                gene_id = tran_id.rsplit(".", 1)[0]

            if tran_id and gene_id:
                # Store multiple variants of the transcript ID
                tran_to_gene[tran_id] = gene_id
                tran_to_gene["rna-" + tran_id] = gene_id
                tran_to_gene["transcript:" + tran_id] = gene_id
                tran_to_gene[tran_id.replace("rna-", "")] = gene_id

    seen = set()
    internal_stop_genes = []

    def flush(gene_id, seq_lines, fout):
        if gene_id in seen:
            print(f"WARNING: Duplicate gene ID skipped: {gene_id}", flush=True)
            return
        # gffread represents the true terminal stop codon as a trailing '*' —
        # expected, and stripped below. A '*' anywhere else means the CDS has a
        # premature stop (bad gene model: assembly/annotation error, or mixed
        # annotation pipelines between species).
        seq = "".join(seq_lines).replace(".", "")
        body = seq[:-1] if seq.endswith("*") else seq
        if "*" in body:
            internal_stop_genes.append(gene_id)
            if INTERNAL_STOP_ACTION == "drop":
                return
            body = body.replace("*", "")
        seen.add(gene_id)
        fout.write(f">{gene_id}\\n")
        fout.write(body + "\\n")

    with open("${fasta}") as fin, open("${meta.id}.clean.fasta", "w") as fout:
        pending = None
        for line in fin:
            if line.startswith(">"):
                if pending is not None:
                    flush(pending[0], pending[1], fout)
                seq_id = line[1:].strip().split()[0]
                gene_id = tran_to_gene.get(seq_id, seq_id)
                # Last resort: strip .tN suffix for Braker style
                if gene_id == seq_id and "." in seq_id:
                    gene_id = seq_id.rsplit(".", 1)[0]
                pending = (gene_id, [])
            elif pending is not None:
                pending[1].append(line.rstrip("\\n"))
        if pending is not None:
            flush(pending[0], pending[1], fout)

    if internal_stop_genes:
        verb = "dropped" if INTERNAL_STOP_ACTION == "drop" else "kept (internal stop codon(s) stripped)"
        print(f"WARNING: {len(internal_stop_genes)} gene(s) had a premature stop codon and were {verb}", flush=True)
        with open("${meta.id}.internal_stop_codons.tsv", "w") as rep:
            rep.write("gene_id\\taction\\n")
            for g in internal_stop_genes:
                rep.write(f"{g}\\t{INTERNAL_STOP_ACTION}\\n")
    EOF
    """
}

