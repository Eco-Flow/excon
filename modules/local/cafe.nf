process CAFE {
    label 'process_low'
    label 'error_ignore'
    tag "cafe"
    container= 'ecoflowucl/cafe:r-4.3.1'

    input:
    path Table
    path tree_newick

    output:
    path("Significant_trees.tre") , emit: cafe_significant_nexus
    path("Out_cafe") , emit: result
    path("Out_cafe_k3") , emit: result_k3
    path("Out_cafe_p_k3") , emit: result_p_k3
    path("N0.tsv") , emit: N0_table
    path "versions.yml", emit: versions

    script:
    """
    cp  ${tree_newick} pruned_tree
    sed -i 's/.prot.fa.largestIsoform//g' pruned_tree 
    sed -i 's/.prot.fa.largestIsoform//g' N0.tsv

    ${projectDir}/bin/cafe_prep.R 

    cafe5 -i hog_gene_counts.tsv -t SpeciesTree_rooted_ultra.txt --cores ${task.cpus} -o Out_cafe 
    cafe5 -i hog_gene_counts.tsv -t SpeciesTree_rooted_ultra.txt --cores ${task.cpus} -o Out_cafe_k3 -k 3
    cafe5 -i hog_gene_counts.tsv -t SpeciesTree_rooted_ultra.txt --cores ${task.cpus} -o Out_cafe_p_k3 -p -k 3

    echo '#nexus\nbegin trees;' > Significant_trees.tre
    grep "*" Out_cafe_k3/Gamma_asr.tre >> Significant_trees.tre
    echo "end;">>Significant_trees.tre

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        R version: \$(R --version | grep "R version" | sed 's/[(].*//' | sed 's/ //g' | sed 's/[^0-9]*//')
    END_VERSIONS
    """
}
