process ORTHOFINDER {
    //tag "$meta.id"
    label 'process_high'
    label 'process_long'
    label 'process_high_memory'
    container = 'biocontainers/orthofinder:3.1.0--hdfd78af_0'

    input:
    tuple val(meta), path(fastas, stageAs: 'input/')

    output:
    tuple val(meta), path("$meta")                                , emit: orthofinder
    path("$meta/Orthogroups/Orthogroups.tsv")                     , emit: orthologues
    path("$meta/Species_Tree/SpeciesTree_rooted_node_labels.txt") , emit: speciestree
    path("$meta/Phylogenetic_Hierarchical_Orthogroups/N1.tsv")    , emit: no_ortho
    path "versions.yml"                                           , emit: versions

    when:
    task.ext.when == null || task.ext.when

    script:
    def args    = task.ext.args ?: ''
    //prefix      = task.ext.prefix ?: "${meta.id}"
    """
    mkdir temp_pickle

    orthofinder \\
        $args \\
        -t $task.cpus \\
        -a $task.cpus \\
        -p temp_pickle \\
        -f input \\
        -n $meta

    mv \\
        input/OrthoFinder/Results_$meta \\
        $meta

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        orthofinder: \$(orthofinder -h | sed -n 's/.*version \\(.*\\) Copy.*/\\1/p')
    END_VERSIONS
    """

    stub:
    def args    = task.ext.args ?: ''
    //prefix      = task.ext.prefix ?: "${meta.id}"
    """
    mkdir -p    $meta/Comparative_Genomics_Statistics
    mkdir       $meta/Gene_Duplication_Events
    mkdir       $meta/Gene_Trees
    mkdir       $meta/Orthogroup_Sequences
    mkdir       $meta/Orthogroups
    mkdir       $meta/Orthologues
    mkdir       $meta/Phylogenetic_Hierarchical_Orthogroups
    mkdir       $meta/Phylogenetically_Misplaced_Genes
    mkdir       $meta/Putative_Xenologs
    mkdir       $meta/Resolved_Gene_Trees
    mkdir       $meta/Single_Copy_Orthologue_Sequences
    mkdir       $meta/Species_Tree
    mkdir       $meta/WorkingDirectory

    touch       $meta/Log.txt

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        orthofinder: \$(orthofinder -h | sed -n 's/.*version \\(.*\\) Copy.*/\\1/p')
    END_VERSIONS
    """
}
