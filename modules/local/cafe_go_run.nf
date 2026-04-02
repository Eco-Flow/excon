process CAFE_GO_RUN {
    tag "${meta.id}"
    label 'process_low'
    container 'ecoflowucl/chopgo:r-4.3.2_python-3.10_perl-5.38'

    input:
    tuple val(meta), path(target), path(background), path(og_go)

    output:
    tuple val(meta), path("*_TopGo_results_ALL.tab"), emit: topgo_results, optional: true
    tuple val(meta), path("*.pdf"),                   emit: pdfs,          optional: true
    tuple val("${task.process}"), val('perl'), eval("perl --version 2>&1 | grep 'version' | sed 's/.*(//; s/[)].*//'"), emit: versions_cafe_perl, topic: versions
    tuple val("${task.process}"), val('R'), eval("R --version 2>&1 | grep 'R version' | sed 's/.*R version //; s/ .*//'"), emit: versions_cafe_R, topic: versions

    script:
    """
    ${projectDir}/bin/cafe_go_run.pl \\
        ${target} \\
        ${background} \\
        OG_GO_format.tsv \\
        ${params.go_cutoff} \\
        ${params.go_type} \\
        ${params.go_max_plot} \\
        ${params.go_algo}
    """
}
