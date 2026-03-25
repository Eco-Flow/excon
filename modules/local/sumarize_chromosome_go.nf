process SUMMARIZE_CHROMO_GO {
    tag "$meta.id"

    container "${ workflow.containerEngine == 'singularity' && !task.ext?.singularity_pull_docker_container ?
        'docker://rocker/tidyverse:4.3.2' :
        'rocker/tidyverse:4.3.2' }"

    input:
    tuple val(meta), path(res_dir)

    output:
    tuple val(meta), path("*.pdf"), emit: plots
    tuple val(meta), path("*.csv"), emit: tables

    script:
    """
    summarize_chromosome_go.R --input ${res_dir}
    """
}
