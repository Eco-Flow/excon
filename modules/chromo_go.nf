process CHROMO_GO {
    label 'chromo_go'
    tag "$sample_id"
    container = 'chriswyatt/chopgo'
    publishDir "$params.outdir/Chromo_Go" , mode: "copy"

    input:

        path gffs
	path goes
	path "Orthogroups.tsv"

    output:

        path( "Unfiltered_Go*" ), emit: chromosome_go_unfilt
	path( "Filtered_dup_Go_*" ), emit: chromosome_go_filt

    script:
    """
    # Run GO chromosome script:
	go_chromosome.pl
    """
}


    
