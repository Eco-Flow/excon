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
    tuple val(meta), path("${meta.id}.go.txt"), emit: go_file

    script:
    """
    python3 <<EOF
    # Build transcript -> gene mapping from GFF
    tran_to_gene = {}
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

            # Get transcript ID - try ID attribute
            tran_id = attrs.get("ID", "")
            # Strip common prefixes
            tran_id_clean = tran_id.replace("rna-", "").replace("transcript:", "")

            # Get gene ID - try Parent attribute
            gene_id = attrs.get("Parent", "")
            if ":" in gene_id:
                gene_id = gene_id.split(":")[-1]
            gene_id = gene_id.replace("gene-", "")

            # Also handle maker/augustus style (gene is in Name or gene attribute)
            if not gene_id:
                gene_id = attrs.get("gene", attrs.get("Name", tran_id_clean))

            # Handle Braker/augustus style where ID is like g1.t1
            # gene ID would be g1, transcript g1.t1
            if "." in tran_id_clean and not gene_id:
                gene_id = tran_id_clean.rsplit(".", 1)[0]

            if tran_id and gene_id:
                tran_to_gene[tran_id] = gene_id
                tran_to_gene[tran_id_clean] = gene_id
                # Also store with rna- prefix variant
                tran_to_gene["rna-" + tran_id_clean] = gene_id
                tran_to_gene["transcript:" + tran_id_clean] = gene_id

    # Parse eggnogmapper annotations and output gene-level GO file
    with open("${annotations}") as f, open("${meta.id}.go.txt", "w") as out:
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

            # Try to get gene ID from mapping, fallback to query itself
            gene_id = tran_to_gene.get(query, None)
            if not gene_id:
                # Try stripping prefixes
                query_clean = query.replace("rna-", "").replace("transcript:", "")
                gene_id = tran_to_gene.get(query_clean, query_clean)
                # Handle g1.t1 style
                if "." in gene_id:
                    gene_id = gene_id.rsplit(".", 1)[0]

            for go in gos.split(","):
                go = go.strip()
                if go.startswith("GO:"):
                    out.write(f"{gene_id}\\t{go}\\n")
    EOF
    """
}
