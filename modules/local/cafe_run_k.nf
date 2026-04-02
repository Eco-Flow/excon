process CAFE_RUN_K {
    tag "k=${k}"
    label 'process_high'
    label 'process_long'
    container 'ecoflowucl/cafe:r-4.3.1' //R version 4.3.1, cafe version 4.2.1 (confusing)

    input:
    path hog_counts
    path species_tree
    path error_model
    each k   // <-- this fans out automatically in DSL2

    output:
    tuple val(k), path("Out_cafe_k${k}/"), emit: results
    path "cafe_k${k}.log",                 emit: log
    tuple val("${task.process}"), val('cafe'), val('4.2.1'), emit: versions_cafe, topic: versions

    script:
    def k_flag = k > 1 ? "-k ${k}" : ""
    def e_flag = error_model.size() > 0 ? "-e${error_model}" : ""  // no space after -e
    """
    cafe5 \\
        -i ${hog_counts} \\
        -t ${species_tree} \\
        --cores ${task.cpus} \\
        ${k_flag} \\
        ${e_flag} \\
        -o Out_cafe_k${k} \\
        2>&1 | tee cafe_k${k}.log
    cafe5_exit=\${PIPESTATUS[0]}
    exit \$cafe5_exit
    """

}
