process CAFE_RUN_BEST {
    tag "${use_poisson ? "k${best_k}_poisson" : "k${best_k}"}"
    label 'process_high'
    label 'process_long'
    container 'ecoflowucl/cafe:r-4.3.1' 

    input:
    path  hog_counts
    path  species_tree
    path  error_model
    val   best_k
    each  use_poisson

    output:
    path "Out_cafe_${use_poisson ? "k${best_k}_poisson" : "k${best_k}"}/", emit: results
    path "cafe_${use_poisson ? "k${best_k}_poisson" : "k${best_k}"}.log",  emit: log
    tuple val("${task.process}"), val('cafe'), val('4.2.1'),                emit: versions_cafe, topic: versions

    script:
    def run_label = use_poisson ? "k${best_k}_poisson" : "k${best_k}"
    def k_flag    = best_k > 1  ? "-k ${best_k}" : ""
    def p_flag    = use_poisson ? "-p" : ""
    def e_flag    = error_model.size() > 0 ? "-e${error_model}" : ""
    """
    cafe5 \\
        -i ${hog_counts} \\
        -t ${species_tree} \\
        --cores ${task.cpus} \\
        ${k_flag} \\
        ${p_flag} \\
        ${e_flag} \\
        -o Out_cafe_${run_label} \\
        2>&1 | tee cafe_${run_label}.log
    cafe5_exit=\${PIPESTATUS[0]}

    exit \$cafe5_exit
    """

    stub:
    def run_label = use_poisson ? "k${best_k}_poisson" : "k${best_k}"
    """
    mkdir -p Out_cafe_${run_label}
    touch Out_cafe_${run_label}/Base_results.txt
    touch cafe_${run_label}.log
    """
}
