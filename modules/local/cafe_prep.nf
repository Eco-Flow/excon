process CAFE_PREP {
    label 'process_high'
    label 'process_long'
    container 'ecoflowucl/cafe:r-4.3.1' //R version 4.3.1 , but cafe version 4.2.1 (confusing)
    
    errorStrategy { 
        if (task.attempt <= 1 && task.exitStatus == 1) {
            log.warn "CAFE_PREP: Base run failed on attempt ${task.attempt} — retrying with differential filtering"
            return 'retry'
        }
        return 'ignore'
    }
    maxRetries 1

    input:
    path table
    path tree_newick

    output:
    path("hog_gene_counts.tsv"),         emit: prepared_counts
    path("SpeciesTree_rooted_ultra.txt"),emit: prepared_tree
    path("N0.tsv"),                      emit: N0_table
    path("Out_cafe"),                    emit: results
    path("Out_cafe/Base_count.tab"),     emit: result_nftest
    path("hog_filtering_report.tsv"),    emit: filtering_report, optional: true
    tuple val("${task.process}"), val('R'), eval("R --version 2>&1 | grep 'R version' | sed 's/R version \\([0-9.]*\\).*/\\1/'"), emit: versions_R, topic: versions
    tuple val("${task.process}"), val('cafe'), val('4.2.1'), emit: versions_cafe, topic: versions

    script:
    def max_differential = params.cafe_max_differential ?: 50
    def use_filtering    = task.attempt > 1
    """
    export PATH=\$PATH:/usr/bin
    set -e

    mv Orthogroups.tsv N0.tsv
    cp ${tree_newick} pruned_tree
    sed -i 's/\\.clean//g' pruned_tree
    sed -i 's/\\.clean//g' N0.tsv

    if [ "${use_filtering}" = "true" ]; then
        echo "CAFE_PREP attempt ${task.attempt}: applying differential filtering (threshold: ${max_differential})"
        /usr/bin/Rscript ${projectDir}/bin/cafe_prep_filtered.R ${max_differential}
    else
        echo "CAFE_PREP attempt ${task.attempt}: no filtering"
        /usr/bin/Rscript ${projectDir}/bin/cafe_prep.R
    fi

    cafe5 \\
        -i hog_gene_counts.tsv \\
        -t SpeciesTree_rooted_ultra.txt \\
        --cores ${task.cpus} \\
        -o Out_cafe \\
        2>&1 | tee cafe_base.log
    cafe5_exit=\${PIPESTATUS[0]}

    exit \$cafe5_exit
    """
}
