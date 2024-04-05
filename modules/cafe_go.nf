process CAFE_GO {
    label 'cafe_go'
    tag "cafe_go"
    publishDir "$params.outdir/Cafe_go/", mode: "copy"
    errorStrategy = 'ignore'    
    container= 'chriswyatt/chopgo'

    input:
        path Table
        path N0_table
        path Go

    output:
        path("*.pdf") , emit: cafe_go_pdfs
	path("CAFE_summary.txt"), emit: cafe_summary

    script:
    """
	cafe_go.pl
	sum_cafe.pl
	plotting_go.R
    """
}
