process CAFE_RUN_K {
    tag "k=${k}"
    label 'process_high'
    label 'process_long'
    container 'ecoflowucl/cafe:r-4.3.1' // R 4.3.1; CAFE5 1.1.0 (cafe5 has no --version: from /opt/CAFE5/config.h)

    input:
    path hog_counts
    path species_tree
    path error_model
    each k

    output:
    tuple val(k), path("Out_cafe_k${k}/"), emit: results
    path "cafe_k${k}.log",                 emit: log
    tuple val("${task.process}"), val('cafe5'), val('1.1.0'), emit: versions_cafe, topic: versions

    script:
    def k_flag = k > 1 ? "-k ${k}" : ""
    def e_flag = error_model.size() > 0 ? "-e${error_model}" : ""  // no space after -e
    // Family-wide significance (--cafe_family_pvalue). Under the gamma model CAFE5
    // writes branch p-values only for families significant at this level.
    def P_flag = "--pvalue ${params.cafe_family_pvalue}"
    // Parameters reused from an earlier run (--cafe_fixed_lambda / --cafe_fixed_alpha)
    def l_flag = params.cafe_fixed_lambda != null ? "--fixed_lambda ${params.cafe_fixed_lambda}" : ""
    def a_flag = (k > 1 && params.cafe_fixed_alpha != null) ? "--fixed_alpha ${params.cafe_fixed_alpha}" : ""
    """
    cafe5 \\
        -i ${hog_counts} \\
        -t ${species_tree} \\
        --cores ${task.cpus} \\
        ${k_flag} \\
        ${e_flag} \\
        ${P_flag} \\
        ${l_flag} \\
        ${a_flag} \\
        -o Out_cafe_k${k} \\
        2>&1 | tee cafe_k${k}.log
    cafe5_exit=\${PIPESTATUS[0]}
    exit \$cafe5_exit
    """

}
