process EGGNOG_TO_GO {
    tag "$meta.id"
    label 'process_single'

    container "${ workflow.containerEngine == 'singularity' && !task.ext.singularity_pull_docker_container ?
        'https://depot.galaxyproject.org/singularity/eggnog-mapper:2.1.13--pyhdfd78af_2' :
        'biocontainers/eggnog-mapper:2.1.13--pyhdfd78af_2' }"

    input:
        tuple val(meta), path(annotations)
        tuple val(meta2), path(gff)

    output:
        tuple val(meta), path("${meta.id}.go.txt"),         emit: go_file
        tuple val(meta), path("${meta.id}.isoform_go.txt"), emit: isoform_go_file
        tuple val("${task.process}"), val('python'), eval("python3 --version | sed 's/Python //'"), emit: versions, topic: versions

    script:
    """
    python3 <<'EOF'
    # Build transcript -> gene AND gene -> [all transcripts] mappings from GFF
    tran_to_gene = {}
    gene_to_tran = {}

    with open("${gff}") as f:
        for line in f:
            if line.startswith("#"):
                continue
            parts = line.strip().split("\\t")
            if len(parts) < 9 or parts[2] != "mRNA":
                continue

            attrs = {}
            for a in parts[8].split(";"):
                if "=" in a:
                    k, v = a.split("=", 1)
                    attrs[k.strip()] = v.strip()

            tran_id = attrs.get("ID", "")
            tran_id_clean = tran_id.replace("rna-", "").replace("transcript:", "")

            gene_id = attrs.get("Parent", "")
            if ":" in gene_id:
                gene_id = gene_id.split(":")[-1]
            gene_id = gene_id.replace("gene-", "")

            if not gene_id:
                gene_id = attrs.get("gene", attrs.get("Name", tran_id_clean))

            if "." in tran_id_clean and not gene_id:
                gene_id = tran_id_clean.rsplit(".", 1)[0]

            if tran_id and gene_id:
                tran_to_gene[tran_id] = gene_id
                tran_to_gene[tran_id_clean] = gene_id
                tran_to_gene["rna-" + tran_id_clean] = gene_id
                tran_to_gene["transcript:" + tran_id_clean] = gene_id

                # Store all isoforms per gene
                gene_to_tran.setdefault(gene_id, []).append(tran_id_clean)

    with open("${annotations}") as f, \
         open("${meta.id}.go.txt", "w") as out_gene, \
         open("${meta.id}.isoform_go.txt", "w") as out_iso:

        for line in f:
            if line.startswith("#"):
                continue

            cols = line.strip().split("\\t")
            if len(cols) < 10:
                continue

            query = cols[0]
            gos = cols[9]

            if gos == "-" or gos == "":
                continue

            gene_id = tran_to_gene.get(query, None)

            if not gene_id:
                query_clean = query.replace("rna-", "").replace("transcript:", "")
                gene_id = tran_to_gene.get(query_clean, query_clean)

                if "." in gene_id:
                    gene_id = gene_id.rsplit(".", 1)[0]

            for go in gos.split(","):
                go = go.strip()
                if go.startswith("GO:"):
                    out_gene.write(f"{gene_id}\\t{go}\\n")

                    # Propagate GO to all isoforms of this gene
                    for tran_id in gene_to_tran.get(gene_id, [gene_id]):
                        out_iso.write(f"{tran_id}\\t{go}\\n")
    EOF
    """
}

