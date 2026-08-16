process DATE_TREE {
    label 'process_single'

    container 'ecoflowucl/cafe:r-4.3.1' //R version 4.3.1, includes ape

    input:
    path tree_newick
    path calibrations

    output:
    path 'SpeciesTree_dated.nwk'   , emit: dated_tree
    path 'dating_calibrations.tsv' , emit: calibrations
    path 'dating_qc.tsv'           , emit: qc
    tuple val("${task.process}"), val('ape'), eval("Rscript -e 'cat(as.character(packageVersion(\"ape\")))'"), emit: versions, topic: versions

    when:
    task.ext.when == null || task.ext.when

    script:
    def model     = params.chronos_model ?: 'discrete'
    def lambda    = params.chronos_lambda ?: 1
    def rate_cats = params.chronos_rate_categories ?: 10
    """
    Rscript ${projectDir}/bin/date_tree.R \\
        ${tree_newick} \\
        ${calibrations} \\
        ${model} \\
        ${lambda} \\
        ${rate_cats}
    """

    stub:
    """
    touch SpeciesTree_dated.nwk
    touch dating_calibrations.tsv
    touch dating_qc.tsv
    """
}
