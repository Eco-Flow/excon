// CAFE process with intelligent retry on differential errors
// First attempt: standard CAFE without filtering
// Second attempt (on failure): CAFE with differential filtering

process CAFE {
    label 'process_medium'
    label 'process_long'
    container 'ecoflowucl/cafe:r-4.3.1' //R version 4.3.1 , but cafe version 4.2.1 (confusing)
    
    errorStrategy { task.attempt <= 1 ? 'retry' : 'ignore' }
    maxRetries 1

    input:
    path Table
    path tree_newick

    output:
    path("Significant_trees.tre") , emit: cafe_significant_nexus, optional: true
    path("Out_cafe") , emit: result, optional: true
    path("Out_cafe_k3") , emit: result_k3, optional: true
    path("Out_cafe_p_k3") , emit: result_p_k3, optional: true
    path("N0.tsv") , emit: N0_table
    path("Out_cafe/Base_count.tab") , emit: result_nftest
    path("hog_filtering_report.tsv"), emit: filtering_report, optional: true
    tuple val("${task.process}"), val('R'), eval ("R --version 2>&1 | grep 'R version' | sed 's/R version \\([0-9.]*\\).*/\\1/'"), emit: versions_R, topic: versions
    tuple val("${task.process}"), val('cafe'), val('4.2.1'), emit: versions_cafe, topic: versions

    script:
    def max_differential = params.cafe_max_differential ?: 50
    def use_filtering = task.attempt > 1 ? true : false
    
    """
    export PATH=\$PATH:/usr/bin
    set -e
    mv Orthogroups.tsv N0.tsv    
    cp ${tree_newick} pruned_tree
    sed -i 's/\\.clean//g' pruned_tree
    sed -i 's/\\.clean//g' N0.tsv

    # Check if this is a retry attempt
    if [ "${use_filtering}" = "true" ]; then
        echo "=================================================="
        echo "CAFE RETRY: Applying differential filtering"
        echo "Threshold: ${max_differential}"
        echo "=================================================="
        
        /usr/bin/Rscript ${projectDir}/bin/cafe_prep_filtered.R ${max_differential}
    else
        echo "=================================================="
        echo "CAFE: First attempt without filtering"
        echo "=================================================="
        
        /usr/bin/Rscript ${projectDir}/bin/cafe_prep.R
    fi

    echo "Running CAFE5 analysis..."
    
    if ! cafe5 -i hog_gene_counts.tsv -t SpeciesTree_rooted_ultra.txt --cores ${task.cpus} -o Out_cafe 2>&1 | tee cafe_base.log; then
        echo "ERROR: CAFE5 base run failed"
        if grep -q "largest size differential" cafe_base.log; then
            echo "Detected size differential error - will retry with filtering"
        fi
        exit 1
    fi
    
    if ! cafe5 -i hog_gene_counts.tsv -t SpeciesTree_rooted_ultra.txt --cores ${task.cpus} -o Out_cafe_k3 -k 3 2>&1 | tee cafe_k3.log; 
then
        echo "ERROR: CAFE5 k3 run failed"
        exit 1
    fi
    
    if ! cafe5 -i hog_gene_counts.tsv -t SpeciesTree_rooted_ultra.txt --cores ${task.cpus} -o Out_cafe_p_k3 -p -k 3 2>&1 | tee cafe_p_k3.log; then
        echo "ERROR: CAFE5 p k3 run failed"
        exit 1
    fi

    echo "CAFE5 completed successfully!"
    
    echo '#nexus' > Significant_trees.tre
    echo 'begin trees;' >> Significant_trees.tre
    grep "*" Out_cafe_k3/Gamma_asr.tre >> Significant_trees.tre || echo "No significant trees found"
    echo "end;" >> Significant_trees.tre

    if [ "${use_filtering}" = "true" ]; then
        echo "CAFE completed successfully on retry with differential filtering (threshold: ${max_differential})"
    else
        echo "CAFE completed successfully on first attempt"
    fi
    """
}
