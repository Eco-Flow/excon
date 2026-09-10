process PLOT_CAFE_GO {
    tag "$tag"
    label 'process_single'
    container "${ workflow.containerEngine == 'singularity' && !task.ext?.singularity_pull_docker_container ?
        'docker://rocker/tidyverse:4.3.2' :
        'rocker/tidyverse:4.3.2' }"

    input:
    tuple val(tag), path(pos_tsv), path(neg_tsv)

    output:
    // ** so these also catch plotting_go_summary.R's outputs, which now live
    // under heatmaps/, dotplots/ and tables/ rather than flat in the task dir.
    tuple val(tag), path("**/*.pdf"),  emit: plots,    optional: true
    tuple val(tag), path("**/*.png"),  emit: pngs,     optional: true
    tuple val(tag), path("**/*.svg"),  emit: svgs,     optional: true
    tuple val(tag), path("**/*.tsv"),  emit: tables,   optional: true
    tuple val(tag), path("*.R"),       emit: rscripts, optional: true
    tuple val(tag), path("README.md"), emit: readme,   optional: true
    tuple val("${task.process}"), val('R'), eval("R --version 2>&1 | grep 'R version' | sed 's/R version \\([0-9.]*\\).*/\\1/'"), emit: versions_R, topic: versions

    script:
    """
    cp ${projectDir}/bin/plotting_go.R .
    cp ${projectDir}/bin/plotting_go_summary.R .
    plotting_go.R
    plotting_go_summary.R
    """
}
