process CAFE_MODEL_COMPARE {
    tag "cafe_model_compare"
    label 'process_low'
    container 'ecoflowucl/cafe:r-4.3.1'    

    input:
    path base_results          // Out_cafe        from CAFE_PREP
    tuple val(meta_g),   path(gamma_results)      // Out_gamma       from CAFE_RUN
    tuple val(meta_pg),  path(gamma_pf_results)   // Out_gamma_pf    from CAFE_RUN

    output:
    path("cafe_model_comparison.tsv"),  emit: comparison_table
    path("Significant_trees.tre"),      emit: significant_trees, optional: true
    path("best_model.txt"),             emit: best_model
    tuple val("${task.process}"), val('R'), eval("R --version 2>&1 | grep 'R version' | sed 's/R version \\([0-9.]*\\).*/\\1/'"), emit: versions_R, topic: versions

    script:
    """
    /usr/bin/Rscript ${projectDir}/bin/cafe_model_compare.R \\
        ${base_results} \\
        ${gamma_results} \\
        ${gamma_pf_results}

    # Build Significant_trees.tre from whichever model won
    best=\$(cat best_model.txt)

    echo '#nexus'          > Significant_trees.tre
    echo 'begin trees;'   >> Significant_trees.tre

    if [ "\$best" = "gamma" ]; then
        grep "*" ${gamma_results}/Gamma_asr.tre    >> Significant_trees.tre || true
    elif [ "\$best" = "gamma_per_family" ]; then
        grep "*" ${gamma_pf_results}/Gamma_asr.tre >> Significant_trees.tre || true
    else
        grep "*" ${base_results}/Base_asr.tre      >> Significant_trees.tre || true
    fi

    echo 'end;' >> Significant_trees.tre
    """
}
