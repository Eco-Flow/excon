process ROOT_TREE {
    tag "$meta.id"
    label 'process_single'

    container "${ workflow.containerEngine == 'singularity' && !task.ext.singularity_pull_docker_container ?
        'https://depot.galaxyproject.org/singularity/eggnog-mapper:2.1.13--pyhdfd78af_2' :
        'biocontainers/eggnog-mapper:2.1.13--pyhdfd78af_2' }"

    input:
    tuple val(meta), path(tree_newick)

    output:
    path 'SpeciesTree_rooted.nwk'              , emit: tree
    tuple val("${task.process}"), val('python'), eval("python3 --version | sed 's/Python //'"), emit: versions, topic: versions

    when:
    task.ext.when == null || task.ext.when

    script:
    def outgroup = params.iqtree_outgroup ? "-g '${params.iqtree_outgroup}'" : ''
    """
    root_tree.py \\
        -i ${tree_newick} \\
        -o SpeciesTree_rooted.nwk \\
        $outgroup
    """

    stub:
    """
    touch SpeciesTree_rooted.nwk
    """
}
