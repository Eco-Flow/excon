process ORTHOFINDER {
    label 'process_medium'
    label 'process_long'
    container = 'ecoflowucl/orthofinder:2.5.5'
    
    input:
    path '*'

    output:
    path("Orthogroups/Orthogroups.tsv") , emit: orthologues
    path("Species_Tree/SpeciesTree_rooted_node_labels.txt") , emit: speciestree
    path("Phylogenetic_Hierarchical_Orthogroups/N0.tsv"), emit: no_ortho
    path "versions.yml", emit: versions

    script:
    """
    if [ -f *.gz ]; then
       gunzip *.gz
    fi

    orthofinder -f . -o My_result
    
    #Remove Dated folder system that prevents unit testing
    mv My_result/*/* .

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        Python version: \$(python --version | cut -f 2 -d " ")
        Orthofinder version: \$(orthofinder | grep version | cut -f 3 -d " ")
    END_VERSIONS
    """
}
