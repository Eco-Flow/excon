process CAFE_RUN_LARGE {
    tag "${hog_counts.baseName}"
    label 'process_single'
    container 'ecoflowucl/cafe:r-4.3.1'

    input:
    path hog_counts     // one family's own split of hog_gene_counts_large.tsv
    path species_tree
    path error_model

    output:
    path "Out_cafe_large_${hog_counts.baseName}/", emit: results, optional: true
    path "cafe_large_${hog_counts.baseName}.log",  emit: log
    tuple val("${task.process}"), val('cafe'), val('4.2.1'), emit: versions_cafe, topic: versions

    script:
    def e_flag = error_model.size() > 0 ? "-e${error_model}" : ""
    def z_flag = params.cafe_zero_root ? "-z" : ""
    def outdir = "Out_cafe_large_${hog_counts.baseName}"
    def logf   = "cafe_large_${hog_counts.baseName}.log"
    """
    # A shared lambda across every high-differential family at once is what made
    # this module never converge (see git history / CAFE5 discussion #132): no
    # single rate can explain both a modest and an extreme size differential in
    # the same fit. Each family now runs on its own, free to find whatever
    # lambda IT needs — a single family imposes no cross-family conflict, so
    # CAFE5 almost always finds a finite-likelihood fit.
    cafe5 \\
        -i ${hog_counts} \\
        -t ${species_tree} \\
        --cores ${task.cpus} \\
        ${e_flag} \\
        ${z_flag} \\
        -o ${outdir} \\
        2>&1 | tee ${logf} || true

    if [ ! -f ${outdir}/Base_results.txt ] || \\
       ! grep -q "Final Likelihood" ${outdir}/Base_results.txt || \\
       grep -q "inf" ${outdir}/Base_results.txt; then
        echo "WARNING: ${hog_counts.baseName} did not converge even on its own — dropping it from the large-family merge." >&2
        rm -rf ${outdir}
    fi
    """

    stub:
    def outdir = "Out_cafe_large_${hog_counts.baseName}"
    """
    mkdir -p ${outdir}
    touch ${outdir}/Base_results.txt
    touch ${outdir}/Base_asr.tre
    touch ${outdir}/Base_branch_probabilities.tab
    touch ${outdir}/Base_change.tab
    touch ${outdir}/Base_count.tab
    touch ${outdir}/Base_family_results.txt
    touch ${outdir}/Base_clade_results.txt
    touch cafe_large_${hog_counts.baseName}.log
    """
}
