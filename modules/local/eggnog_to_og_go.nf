process EGGNOG_TO_OG_GO {
    tag "og_go"
    label 'process_single'

    container "${ workflow.containerEngine == 'singularity' && !(task.ext?.singularity_pull_docker_container) ?
        'https://depot.galaxyproject.org/singularity/eggnog-mapper:2.1.13--pyhdfd78af_2' :
        'biocontainers/eggnog-mapper:2.1.13--pyhdfd78af_2' }"

    input:
    path go_files      // per-species *.go.txt files
    path orthogroups   // Orthogroups.tsv from OrthoFinder

    output:
    path "OG_GO_format.tsv", emit: og_go
    tuple val("${task.process}"), val('python'), eval("python3 --version | sed 's/Python //'"), emit: versions, topic: versions

    script:
    """
    python3 <<EOF
    import os, glob

    # Build gene -> GO mapping from all species go files
    gene_go = {}
    for f in glob.glob("*.go.txt"):
        with open(f) as fh:
            for line in fh:
                parts = line.strip().split("\\t")
                if len(parts) == 2:
                    gene, go = parts
                    gene_go.setdefault(gene, set()).add(go)

    # Build OG -> genes mapping from Orthogroups.tsv
    og_go = {}
    with open("${orthogroups}") as fh:
        header = fh.readline().strip().split("\\t")
        for line in fh:
            parts = line.strip().split("\\t")
            og = parts[0]
            for genes_str in parts[1:]:
                for gene in genes_str.split(", "):
                    gene = gene.strip()
                    if gene and gene in gene_go:
                        og_go.setdefault(og, set()).update(gene_go[gene])

    # Write OG_GO_format.tsv
    with open("OG_GO_format.tsv", "w") as out:
        for og, gos in og_go.items():
            for go in sorted(gos):
                out.write(f"{og}\\t{go}\\n")
    EOF
    """
}
