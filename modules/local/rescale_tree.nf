process RESCALE_TREE {
    label 'process_single'

    container 'ecoflowucl/chopgo:r-4.3.2_python-3.10_perl-5.38'

    input:
    path tree_newick

    output:
    path 'SpeciesTree_rescaled.nwk', emit: rescaled_tree

    script:
    def scale_factor = params.tree_scale_factor ?: 1000
    """
    rescale_tree.py \\
        -i ${tree_newick} \\
        -o SpeciesTree_rescaled.nwk \\
        -s ${scale_factor}

    """
}
