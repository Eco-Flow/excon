process ORTHOFINDER {
    label 'orthofinder'
    publishDir "$params.outdir/Orthofinder/", mode: "copy"
    container = 'davidemms/orthofinder:2.5.4'
    
    input:
        path '*'
               
    output:
        path("My_result/*/Orthogroups/Orthogroups.tsv") , emit: orthologues
	path("My_result/*/Species_Tree/SpeciesTree_rooted_node_labels.txt") , emit:speciestree
	path("My_result/*/Phylogenetic_Hierarchical_Orthogroups/N0.tsv"), emit: no_ortho

    script:
    """
	ulimit -Sn 4096

        count=`ls -1 *.gz 2>/dev/null | wc -l`
        if [ \$count != 0 ]
        then
	    gunzip *.gz
        fi

	    orthofinder -f . -o My_result
    """
}
