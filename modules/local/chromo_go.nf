process CHROMO_GO {
    label 'process_single'
    tag "chromo_go"
    container 'ecoflowucl/chopgo:r-4.3.2_python-3.10_perl-5.38'

    input:
    tuple val(meta), path(gff), path(go)
    path Orthogroups

    output:
    tuple val(meta), path("Filtered_dup_Go_*"), emit: chromosome_go_filt
    tuple val(meta), path("Unfiltered_Go*"), emit: chromosome_go_unfilt
    path "versions.yml", emit: versions

    script:
    """
    # Make R available to perl backtick calls
    ln -s /usr/bin/R ./R
    export PATH=\$PWD:\$PATH
    
    go_chromosome.pl

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        Perl version: \$(perl --version | grep "version" | sed 's/.*(//g' | sed 's/[)].*//')
    END_VERSIONS
    """
}
