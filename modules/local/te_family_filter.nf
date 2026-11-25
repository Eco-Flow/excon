process TE_FAMILY_FILTER {
    tag "te_family_filter"
    label 'process_low'

    container "${ workflow.containerEngine == 'singularity' && !task.ext.singularity_pull_docker_container ?
        'https://depot.galaxyproject.org/singularity/eggnog-mapper:2.1.13--pyhdfd78af_2' :
        'biocontainers/eggnog-mapper:2.1.13--pyhdfd78af_2' }"

    input:
    path annotations   // *.emapper.annotations, one per species
    path orthogroups   // Orthogroups.tsv / N0.tsv
    path signatures    // --cafe_te_signatures

    output:
    path "filtered/*",               emit: orthogroups
    path "te_families.tsv",          emit: report
    path "te_filter_summary.txt",    emit: summary
    tuple val("${task.process}"), val('python'), eval("python3 --version | sed 's/Python //'"), emit: versions, topic: versions

    script:
    """
    python3 ${projectDir}/bin/te_family_filter.py \\
        --orthogroups ${orthogroups} \\
        --annotations ${annotations} \\
        --signatures ${signatures} \\
        --min-fraction ${params.cafe_te_min_fraction} \\
        --out-dir filtered \\
        > te_filter_summary.txt
    cat te_filter_summary.txt
    """
}
