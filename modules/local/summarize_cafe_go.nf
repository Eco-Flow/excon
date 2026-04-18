process SUMMARIZE_CAFE_GO {
    tag "$tag"
    label 'process_single'
    container 'ecoflowucl/chopgo:r-4.3.2_python-3.10_perl-5.38'

    input:
    tuple val(tag), path(topgo_files)

    output:
    tuple val(tag), path("*.pdf"), emit: plots,  optional: true
    tuple val(tag), path("*.tsv"), emit: tables, optional: true
    tuple val("${task.process}"), val('perl'), eval("perl --version 2>&1 | grep 'version' | sed 's/.*(//; s/[)].*//'"), emit: versions_perl, topic: versions
    tuple val("${task.process}"), val('R'),    eval("R --version 2>&1 | grep 'R version' | sed 's/R version \\([0-9.]*\\).*/\\1/'"), emit: versions_R,    topic: versions

    script:
    """
    sum_cafe.pl
    plotting_go.R
    """
}
