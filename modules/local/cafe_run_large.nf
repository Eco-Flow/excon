process CAFE_RUN_LARGE {
    tag "${hog_counts.baseName}"
    label 'process_low'
    container 'ecoflowucl/cafe:r-4.3.1'

    // A handful of families need many more Nelder-Mead iterations than most to
    // pin down their own lambda, and can run past process_low's 2h/4GB (killed
    // externally by the scheduler rather than exiting cleanly, so the script's
    // own non-convergence handling below never gets a chance to run). This track
    // is best-effort — retry with more time/memory a couple of times (task.attempt
    // scales process_low's 2h/4GB up each retry), then drop the family from the
    // merge rather than fail the whole run.
    errorStrategy { task.attempt <= 3 ? 'retry' : 'ignore' }
    maxRetries 3

    input:
    // Bundled into one tuple (built with .combine() in main.nf) rather than three
    // separate positional channels, so species_tree/error_model are unambiguously
    // paired with every one of the 158 hog_counts items rather than relying on
    // Nextflow's implicit broadcast of a singleton channel alongside a multi-item one.
    tuple path(hog_counts), path(species_tree), path(error_model)

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
