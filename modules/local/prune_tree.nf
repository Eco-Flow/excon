process PRUNE_TREE {
    label 'process_single'

    container 'ecoflowucl/cafe:r-4.3.1' //R version 4.3.1, includes ape

    input:
    path tree_newick
    // A file of tip names for --cafe_species, staged so it is visible inside the
    // container's work directory; assets/NO_FILE when the selection is inline instead.
    path species_file

    output:
    path 'SpeciesTree_pruned.nwk'    , emit: tree
    path 'pruned_tree_species.tsv'   , emit: species
    tuple val("${task.process}"), val('ape'), eval("Rscript -e 'cat(as.character(packageVersion(\"ape\")))'"), emit: versions, topic: versions

    when:
    task.ext.when == null || task.ext.when

    script:
    def mode      = params.cafe_clade ? 'clade' : 'species'
    def selection = (species_file.name != 'NO_FILE') ? species_file : (params.cafe_clade ?: params.cafe_species)
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
