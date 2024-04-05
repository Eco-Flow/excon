process GO_EXPANSION {
    label 'go_expansion'
    tag "GO_expansion"
    publishDir "$params.outdir/Go_expansion/", mode: "copy"
    errorStrategy = 'ignore'    
    container = 'chriswyatt/chopgo'

    input:
        path Go_counts

    output:
        path("GO_table_counts.tsv") , emit: go_count_all
	path("GO_table_counts_forCAFE.tsv") , emit: go_count_table

    script:
    """
	#ls -lash ${Go_counts}
	Expansion_summary.pl
    """
}
