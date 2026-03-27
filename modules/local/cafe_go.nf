process CAFE_GO {
    label 'process_low'
    label 'process_high_memory'
    container 'ecoflowucl/chopgo:r-4.3.2_python-3.10_perl-5.38'

    input:
    path Table
    path N0_table
    path Go

    output:
    path("*.pdf") , emit: cafe_go_pdfs
    path("CAFE_summary.txt"), emit: cafe_summary
    tuple val("${task.process}"), val('perl'), eval("perl --version 2>&1 | grep 'version' | sed 's/.*(//; s/[)].*//'"), emit: versions_cafe_go, topic: versions
    tuple val("${task.process}"), val('R'), eval("R --version 2>&1 | grep 'R version' | sed 's/.*R version //; s/ .*//'"),, emit: versions_cafe_go, topic: versions

    script:
    """
    ${projectDir}/bin/cafe_go.pl ${params.go_cutoff} ${params.go_type} ${params.go_max_plot}
    ${projectDir}/bin/sum_cafe.pl
    ${projectDir}/bin/plotting_go.R
    """
}
