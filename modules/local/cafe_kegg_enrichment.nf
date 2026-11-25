process CAFE_KEGG_ENRICHMENT {
    tag "cafe_kegg"
    label 'process_low'

    container "${ workflow.containerEngine == 'singularity' && !task.ext.singularity_pull_docker_container ?
        'https://depot.galaxyproject.org/singularity/eggnog-mapper:2.1.13--pyhdfd78af_2' :
        'biocontainers/eggnog-mapper:2.1.13--pyhdfd78af_2' }"

    input:
    path annotations                         // *.emapper.annotations, one per species
    path n0_table                            // N0.tsv as CAFE_PREP used it
    path cafe_results                        // the selected main model's CAFE5 results
    path sig_changes                         // significant_changes_per_node.tsv
    path pathway_names, stageAs: 'names/*'   // --cafe_kegg_pathway_names, or NO_FILE

    output:
    path "kegg_enrichment.tsv",             emit: table
    path "kegg_enrichment_significant.tsv", emit: significant
    path "kegg_enrichment_summary.txt",     emit: summary
    tuple val("${task.process}"), val('python'), eval("python3 --version | sed 's/Python //'"), emit: versions, topic: versions

    script:
    // Names from the given file, or downloaded from KEGG (skipped, with a warning, offline)
    def names = params.cafe_kegg_pathway_names ? "--pathway-names ${pathway_names}" : '--download-names'
    """
    python3 ${projectDir}/bin/cafe_kegg_enrichment.py \\
        --annotations ${annotations} \\
        --orthogroups ${n0_table} \\
        --cafe-results ${cafe_results} \\
        --sig-changes ${sig_changes} \\
        ${names} \\
        --padj-cutoff ${params.cafe_kegg_padj} \\
        -o kegg_enrichment.tsv \\
        > kegg_enrichment_summary.txt
    cat kegg_enrichment_summary.txt
    """
}
