process GO_EXPANSION {
    label 'process_low'
    label 'error_ignore'
    tag "GO_expansion"
    container = 'ecoflowucl/chopgo:r-4.3.2_python-3.10_perl-5.38'

    input:
    path Go_counts

    output:
    path("GO_table_counts.tsv") , emit: go_count_all
    path("GO_table_counts_forCAFE.tsv") , emit: go_count_table
    path "versions.yml", emit: versions

    script:
    """
    Expansion_summary.pl

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        Perl version: \$(perl --version | grep "version" | sed 's/.*(//g' | sed 's/[)].*//')
    END_VERSIONS
    """
}
