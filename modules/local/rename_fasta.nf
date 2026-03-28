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
    tuple val("${task.process}"), val('python'), eval("python3 --version | sed 's/Python //'"), emit: versions, topic: versions

    script:
    """
    python3 <<EOF
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

    with open("${fasta}") as fin, open("${meta.id}.clean.fasta", "w") as fout:
        for line in fin:
            if line.startswith(">"):
                seq_id = line[1:].strip().split()[0]
                gene_id = tran_to_gene.get(seq_id, seq_id)
                # Last resort: strip .tN suffix for Braker style
                if gene_id == seq_id and "." in seq_id:
                    gene_id = seq_id.rsplit(".", 1)[0]
                fout.write(f">{gene_id}\\n")
            else:
                # Also strip stop codons
                fout.write(line.replace(".", "").replace("*", ""))
    EOF
    """
}
