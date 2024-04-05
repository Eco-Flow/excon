process CAFE {
    label 'cafe'
    tag "cafe"
    publishDir "$params.outdir/Cafe/", mode: "copy"
    errorStrategy = 'ignore'    
    container= 'chriswyatt/cafe_r'

    input:
        path Table
        path tree_newick

    output:
        path("Significant_trees.tre") , emit: cafe_significant_nexus
        path("Out_cafe") , emit: result
        path("Out_cafe_k3") , emit: result_k3
        path("Out_cafe_p_k3") , emit: result_p_k3
	path("N0.ex.tsv") , emit: N0_table

    script:
    """
	cp  ${tree_newick} pruned_tree
	sed -i 's/.prot.fa.largestIsoform//g' pruned_tree 	
	sed -i 's/.prot.fa.largestIsoform//g' N0.tsv
	perl -pe 's/\\r\\n|\\n|\\r/\\n/g' N0.tsv > N0.ex.tsv

	cafe_prep.R 

	cafe5 -i hog_gene_counts.tsv -t SpeciesTree_rooted_ultra.txt --cores 8 -o Out_cafe 
	cafe5 -i hog_gene_counts.tsv -t SpeciesTree_rooted_ultra.txt --cores 8 -o Out_cafe_k3 -k 3
	cafe5 -i hog_gene_counts.tsv -t SpeciesTree_rooted_ultra.txt --cores 8 -o Out_cafe_p_k3 -p -k 3


	echo '#nexus\nbegin trees;' > Significant_trees.tre
	grep "*" Out_cafe_k3/Gamma_asr.tre >> Significant_trees.tre
	echo "end;">>Significant_trees.tre

    """
}
