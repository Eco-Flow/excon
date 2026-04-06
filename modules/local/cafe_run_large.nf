process CAFE_RUN_LARGE {
    label 'process_high'
    label 'process_long'
    container 'ecoflowucl/cafe:r-4.3.1'

    input:
    path  hog_counts_large
    path  species_tree
    path  error_model
    val   lambda

    output:
    path "Out_cafe_large/",   emit: results
    path "cafe_large.log",    emit: log
    tuple val("${task.process}"), val('cafe'), val('4.2.1'), emit: versions_cafe, topic: versions

    script:
    def e_flag = error_model.size() > 0 ? "-e${error_model}" : ""
    """
    # Large-differential families often fail with the estimated lambda.
    # Retry with progressively smaller lambda values as recommended in
    # https://github.com/hahnlab/CAFE5/discussions/132
    converged=false
    for lambda_try in ${lambda} 0.0001 0.00001 0.000001 0.0000001; do
        echo "Trying CAFE_RUN_LARGE with lambda=\${lambda_try}"
        rm -rf Out_cafe_large
        cafe5 \\
            -i ${hog_counts_large} \\
            -t ${species_tree} \\
            --cores ${task.cpus} \\
            -l \${lambda_try} \\
            ${e_flag} \\
            -o Out_cafe_large \\
            2>&1 | tee cafe_large.log || true

        # Check if CAFE5 produced a usable result (finite likelihood)
        if [ -f Out_cafe_large/Base_results.txt ] && \
           grep -q "Final Likelihood" Out_cafe_large/Base_results.txt && \
           ! grep -q "inf" Out_cafe_large/Base_results.txt; then
            echo "Converged with lambda=\${lambda_try}"
            converged=true
            break
        fi
        echo "Did not converge with lambda=\${lambda_try}, trying next..."
    done

    if [ "\${converged}" = "false" ]; then
        echo "WARNING: CAFE_RUN_LARGE did not converge with any lambda — results will be incomplete."
    fi

    # Ensure output directory exists so downstream steps don't fail
    mkdir -p Out_cafe_large
    """

    stub:
    """
    mkdir -p Out_cafe_large
    touch Out_cafe_large/Base_results.txt
    touch cafe_large.log
    """
}
