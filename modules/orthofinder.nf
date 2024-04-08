process ORTHOFINDER {
    label 'process_medium'
    label 'process_long'
    container = 'ecoflowucl/orthofinder:2.5.5'
    
    input:
    path '*'

    output:
    path("My_result/*/Orthogroups/Orthogroups.tsv") , emit: orthologues
    path("My_result/*/Species_Tree/SpeciesTree_rooted_node_labels.txt") , emit:speciestree
    path("My_result/*/Phylogenetic_Hierarchical_Orthogroups/N0.tsv"), emit: no_ortho

    script:
    """
    if [ -f *.gz ]; then
       gunzip *.gz
    fi

    orthofinder -f . -o My_result
    """
}
