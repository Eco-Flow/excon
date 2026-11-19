process CAFE_RUN_ZERO_ROOT {
    tag "cafe_zero_root"
    label 'process_medium'
    container 'ecoflowucl/cafe:r-4.3.1'

    // --cafe_zero_root: a sensitivity analysis of the families the main model dropped
    // for having no genes on one side of the root. Run separately with CAFE5's -z, and
    // never merged into the main results: CAFE5 doesn't recommend -z in general, as
    // these families' root sizes (and so their changes near the root) are poorly
    // determined. Lambda is fixed to the main families' estimate, as for large families,
    // so the results are on the main model's rate scale.

    input:
    path prepared_counts   // hog_gene_counts.tsv: the families given to the main model
    path not_modelled      // families_not_modelled.tsv from CAFE_MODEL_COMPARE
    path species_tree
    path error_model
    val  fixed_lambda

    output:
    path "Out_cafe_zero_root/",        emit: results, optional: true
    path "zero_root_counts.tsv",       emit: counts
    path "cafe_zero_root.log",         emit: log
    tuple val("${task.process}"), val('cafe5'), val('1.1.0'), emit: versions_cafe, topic: versions

    script:
    def e_flag = error_model.size() > 0 ? "-e${error_model}" : ""
    """
    # The dropped families' rows of the main count table
    awk -F'\\t' 'NR == FNR { if (FNR > 1) keep[\$1] = 1; next } FNR == 1 || (\$2 in keep)' \\
        ${not_modelled} ${prepared_counts} > zero_root_counts.tsv

    cafe5 \\
        -i zero_root_counts.tsv \\
        -t ${species_tree} \\
        --cores ${task.cpus} \\
        -z \\
        --fixed_lambda ${fixed_lambda} \\
        ${e_flag} \\
        --pvalue ${params.cafe_family_pvalue} \\
        -o Out_cafe_zero_root \\
        2>&1 | tee cafe_zero_root.log || true

    if [ ! -f Out_cafe_zero_root/Base_results.txt ] || grep -q "inf" Out_cafe_zero_root/Base_results.txt; then
        echo "WARNING: the zero-root sensitivity run has no finite likelihood; see cafe/zero_root/cafe_zero_root.log" >&2
        rm -rf Out_cafe_zero_root
    fi
    """

    stub:
    """
    mkdir -p Out_cafe_zero_root
    touch Out_cafe_zero_root/Base_results.txt zero_root_counts.tsv cafe_zero_root.log
    """
}
