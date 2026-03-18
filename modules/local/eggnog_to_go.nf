process EGGNOG_TO_GO {
    tag "$meta.id"
    label 'process_single'

    conda "conda-forge::python=3.10"
    container "${ workflow.containerEngine == 'singularity' && !task.ext.singularity_pull_docker_container ?
        'https://depot.galaxyproject.org/singularity/python:3.10' :
        'biocontainers/python:3.10' }"

    input:
    tuple val(meta), path(annotations)

    output:
    tuple val(meta), path("${meta.id}.go.txt"), emit: go_file

    script:
    """
    python3 <<EOF
    import sys

    with open("${annotations}") as f, open("${meta.id}.go.txt", "w") as out:
        for line in f:
            if line.startswith("#"):
                continue
            cols = line.strip().split("\\t")
            if len(cols) < 10:
                continue
            gene = cols[0]
            gos  = cols[9]
            if gos == "-" or gos == "":
                continue
            for go in gos.split(","):
                go = go.strip()
                if go.startswith("GO:"):
                    out.write(f"{gene}\\t{go}\\n")
    EOF
    """
}