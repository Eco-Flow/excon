process CHROMO_GO {
    label 'process_single'
    tag "${meta.id}"
    container 'ecoflowucl/chopgo:r-4.3.2_python-3.10_perl-5.38'

    input:
    tuple val(meta), path(gff), path(go)
    path Orthogroups

    output:
    tuple val(meta), path("Filtered_dup_Go_*"), emit: chromosome_go_filt
    tuple val(meta), path("Unfiltered_Go*"), emit: chromosome_go_unfilt
    tuple val("${task.process}"), val('perl'), eval("perl --version 2>&1 | grep 'version' | sed 's/.*(//; s/[)].*//'"), emit: versions_cafe_perl, topic: versions
    tuple val("${task.process}"), val('R'), eval("R --version 2>&1 | grep 'R version' | sed 's/.*R version //; s/ .*//'"), emit: versions_cafe_R, topic: versions

    script:
    """
    # Make R available to perl backtick calls
    ln -s /usr/bin/R ./R
    export PATH=\$PWD:\$PATH

    # go_chromosome.pl expects the orthogroups file to be named Orthogroups.tsv
    [ "${Orthogroups}" = "Orthogroups.tsv" ] || ln -s ${Orthogroups} Orthogroups.tsv

    go_chromosome.pl ${params.go_algo}

    """
}
