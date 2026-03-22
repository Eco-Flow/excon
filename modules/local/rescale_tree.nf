process RESCALE_TREE {
    label 'process_single'

    container "${ workflow.containerEngine == 'singularity' && !task.ext.singularity_pull_docker_container ?
        'https://depot.galaxyproject.org/singularity/eggnog-mapper:2.1.13--pyhdfd78af_2' :
        'biocontainers/eggnog-mapper:2.1.13--pyhdfd78af_2' }"

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
