process PRUNE_TREE {
    label 'process_single'

    container 'ecoflowucl/cafe:r-4.3.1' //R version 4.3.1, includes ape

    input:
    path tree_newick

    output:
    path 'SpeciesTree_pruned.nwk'    , emit: tree
    path 'pruned_tree_species.tsv'   , emit: species
    tuple val("${task.process}"), val('ape'), eval("Rscript -e 'cat(as.character(packageVersion(\"ape\")))'"), emit: versions, topic: versions

    when:
    task.ext.when == null || task.ext.when

    script:
    def mode      = params.cafe_clade ? 'clade' : 'species'
    def selection = params.cafe_clade ?: params.cafe_species
    """
    Rscript ${projectDir}/bin/prune_tree.R \\
        ${tree_newick} \\
        ${mode} \\
        '${selection}'
    """

    stub:
    """
    touch SpeciesTree_pruned.nwk
    touch pruned_tree_species.tsv
    """
}
