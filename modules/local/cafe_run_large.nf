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
    cafe5 \\
        -i ${hog_counts_large} \\
        -t ${species_tree} \\
        --cores ${task.cpus} \\
        -l ${lambda} \\
        ${e_flag} \\
        -o Out_cafe_large \\
        2>&1 | tee cafe_large.log
    cafe5_exit=\${PIPESTATUS[0]}
    exit \$cafe5_exit
    """

    stub:
    """
    mkdir -p Out_cafe_large
    touch Out_cafe_large/Base_results.txt
    touch cafe_large.log
    """
}
