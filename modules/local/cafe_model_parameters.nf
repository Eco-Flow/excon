process CAFE_MODEL_PARAMETERS {
    tag "cafe_model_parameters"
    label 'process_single'
    container 'ecoflowucl/chopgo:r-4.3.2_python-3.10_perl-5.38'

    input:
    path k_results       // Out_cafe_k<k>/ from every CAFE_RUN_K task (uniform root)
    path poisson_results // Out_cafe_k<best>_poisson/ from CAFE_RUN_BEST
    path poisson_log     // its log: CAFE5 writes the fitted Poisson mean only there
    path error_model     // Base_error_model.txt from CAFE_PREP (empty if none)
    path best_model      // best_model.txt from CAFE_MODEL_COMPARE
    path best_k          // best_k.txt from CAFE_SELECT_K

    output:
    path "cafe_model_parameters.tsv", emit: table
    tuple val("${task.process}"), val('python'), val('3.10.13'), emit: versions_python, topic: versions

    script:
    def fixed = [ params.cafe_fixed_lambda != null ? 'lambda'       : null,
                  params.cafe_fixed_alpha  != null ? 'alpha'        : null,
                  params.cafe_error_model          ? 'error_model'  : null,
                  params.cafe_poisson      != null ? 'poisson_mean' : null ].findAll { it }.join(',')
    """
    /usr/bin/python3.10 ${projectDir}/bin/cafe_model_parameters.py \\
        --results ${k_results} ${poisson_results} \\
        --poisson-log ${poisson_log} \\
        --error-model ${error_model} \\
        --best-model ${best_model} \\
        --best-k ${best_k} \\
        --fixed "${fixed}" \\
        -o cafe_model_parameters.tsv
    """

    stub:
    """
    printf "directory\\tk\\troot_distribution\\tpoisson_mean\\tlambda\\talpha\\tepsilon\\tneg_lnL\\tmax_possible_lambda\\tselected\\tfixed_parameters\\n" > cafe_model_parameters.tsv
    """
}
