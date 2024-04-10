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
    path "versions.yml", emit: versions

    script:
    """
    ${projectDir}/bin/cafe_go.pl
    ${projectDir}/bin/sum_cafe.pl
    ${projectDir}/bin/plotting_go.R

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        R version: \$(R --version | grep "R version" | sed 's/[(].*//' | sed 's/ //g' | sed 's/[^0-9]*//')
        Perl version: \$(perl --version | grep "version" | sed 's/.*(//g' | sed 's/[)].*//')
    END_VERSIONS
    """
}
