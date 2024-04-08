process CAFE_GO {
    label 'process_low'
    label 'error_ignore'
    tag "cafe_go"
    container = 'ecoflowucl/chopgo:r-4.3.2_python-3.10_perl-5.38'

    input:
    path Table
    path N0_table
    path Go

    output:
    path("*.pdf") , emit: cafe_go_pdfs
    path("CAFE_summary.txt"), emit: cafe_summary

    script:
    """
    ${projectDir}/bin/cafe_go.pl
    ${projectDir}/bin/sum_cafe.pl
    ${projectDir}/bin/plotting_go.R
    """
}
