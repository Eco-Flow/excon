process SUMMARIZE_CAFE_GO {
    tag "$tag"
    label 'process_single'
    container 'ecoflowucl/chopgo:r-4.3.2_python-3.10_perl-5.38'

    input:
    tuple val(tag), path(topgo_files)

    output:
    tuple val(tag), path("Go_summary_pos.tsv"),          emit: pos_tsv
    tuple val(tag), path("Go_summary_neg.tsv"),          emit: neg_tsv
    tuple val(tag), path("Go_summary_posneg_merged.tsv"), emit: merged_tsv
    tuple val("${task.process}"), val('perl'), eval("perl --version 2>&1 | grep 'version' | sed 's/.*(//; s/[)].*//'"), emit: versions_perl, topic: versions

    script:
    """
    sum_cafe.pl
    """
}

process PLOT_CAFE_GO {
    tag "$tag"
    label 'process_single'
    container "${ workflow.containerEngine == 'singularity' && !task.ext?.singularity_pull_docker_container ?
        'docker://rocker/tidyverse:4.3.2' :
        'rocker/tidyverse:4.3.2' }"

    input:
    tuple val(tag), path(pos_tsv), path(neg_tsv)

    output:
    tuple val(tag), path("*.pdf"), emit: plots, optional: true
    tuple val("${task.process}"), val('R'), eval("R --version 2>&1 | grep 'R version' | sed 's/R version \\([0-9.]*\\).*/\\1/'"), emit: versions_R, topic: versions

    script:
    """
    plotting_go.R
    """
}
