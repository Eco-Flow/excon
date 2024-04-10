process CHROMO_GO {
    label 'process_low'
    tag "$sample_id"
    container = 'ecoflowucl/chopgo:r-4.3.2_python-3.10_perl-5.38'

    input:
    path gffs
    path goes
    path "Orthogroups.tsv"

    output:
    path( "Unfiltered_Go*" ), emit: chromosome_go_unfilt
    path( "Filtered_dup_Go_*" ), emit: chromosome_go_filt
    path "versions.yml", emit: versions

    script:
    """
    go_chromosome.pl

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        Perl version: \$(perl --version | grep "version" | sed 's/.*(//g' | sed 's/[)].*//')
    END_VERSIONS
    """
}
