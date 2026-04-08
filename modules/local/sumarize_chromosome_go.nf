process SUMMARIZE_CHROMO_GO {
    tag "$meta.id"
    label 'process_single'
    container "${ workflow.containerEngine == 'singularity' && !task.ext?.singularity_pull_docker_container ?
        'docker://rocker/tidyverse:4.3.2' :
        'rocker/tidyverse:4.3.2' }"

    input:
    tuple val(meta), path(res_dir)

    output:
    tuple val(meta), path("*.pdf"), emit: plots,  optional: true
    tuple val(meta), path("*.csv"), emit: tables, optional: true
    tuple val("${task.process}"), val('R'), eval ("R --version 2>&1 | grep 'R version' | sed 's/R version \\([0-9.]*\\).*/\\1/'"), emit: versions_R, topic: versions

    script:
    """
    n_files=\$(find -L ${res_dir} -name "*_res.tab" 2>/dev/null | wc -l)
    if [ "\$n_files" -eq 0 ]; then
        echo "No *_res.tab files found in ${res_dir} — skipping summarization (no significant GO terms)."
    else
        summarize_chromosome_go.R --input ${res_dir}
    fi
    """
}
