process SUMMARIZE_CHROMO_GO {
    tag "$meta.id"
    label 'process_single'
    container "${ workflow.containerEngine == 'singularity' && !task.ext?.singularity_pull_docker_container ?
        'docker://rocker/tidyverse:4.3.2' :
        'rocker/tidyverse:4.3.2' }"

    input:
    tuple val(meta), path(res_dir)

    output:
    tuple val(meta), path("*.pdf"), emit: plots
    tuple val(meta), path("*.csv"), emit: tables
    tuple val("${task.process}"), val('R'), eval ("R --version 2>&1 | grep 'R version' | sed 's/R version \\([0-9.]*\\).*/\\1/'"), emit: versions_R, topic: versions    

    script:
    """
    summarize_chromosome_go.R --input ${res_dir}
    """
}
