process MAKE_ULTRAMETRIC {
    label 'process_single'

    container "${ workflow.containerEngine == 'singularity' && !task.ext.singularity_pull_docker_container ?
        'https://depot.galaxyproject.org/singularity/eggnog-mapper:2.1.13--pyhdfd78af_2' :
        'biocontainers/eggnog-mapper:2.1.13--pyhdfd78af_2' }"

    input:
    path tree_newick

    output:
    path 'SpeciesTree_ultrametric.nwk', emit: rescaled_tree
    tuple val("${task.process}"), val('python'), eval("python3 --version | sed 's/Python //'"), emit: versions, topic: versions

    script:
    def root_age = params.tree_scale_factor ?: 1
    """
    make_ultrametric.py \\
        ${tree_newick} \\
        ${root_age} \\
        > SpeciesTree_ultrametric.nwk
    """
}
