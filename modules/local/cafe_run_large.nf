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
    // paired with every hog_counts item rather than relying on Nextflow's implicit
    // broadcast of a singleton channel alongside a multi-item one.
    // fixed_lambda: the lambda to fix (--cafe_large_mode fixed), or '' to estimate
    // the family's own (--cafe_large_mode per_family).
    tuple path(hog_counts), path(species_tree), path(error_model), val(fixed_lambda)

    output:
    path "Out_cafe_large_${hog_counts.baseName}/", emit: results, optional: true
    path "large_status_${hog_counts.baseName}.tsv", emit: status
    path "cafe_large_${hog_counts.baseName}.log",  emit: log
    tuple val("${task.process}"), val('cafe5'), val('1.1.0'), emit: versions_cafe, topic: versions

    script:
    def e_flag = error_model.size() > 0 ? "-e${error_model}" : ""
    def z_flag = params.cafe_zero_root ? "-z" : ""
    // Family-wide significance (--cafe_family_pvalue). Under the gamma model CAFE5
    // writes branch p-values only for families significant at this level.
    def P_flag = "--pvalue ${params.cafe_family_pvalue}"
    def l_flag = fixed_lambda ? "--fixed_lambda ${fixed_lambda}" : ""
    def outdir = "Out_cafe_large_${hog_counts.baseName}"
    def logf   = "cafe_large_${hog_counts.baseName}.log"
    """
    # --cafe_large_mode fixed (default): lambda is fixed to the value estimated on the
    # families of the main run, as the CAFE5 tutorial does for families too large to
    # estimate it from. With lambda fixed, each family's likelihood, reconstruction and
    # p-value don't depend on the others, so running them one at a time gives the same
    # results as one CAFE5 call, and one family that can't be fitted doesn't stop the rest.
    #
    # --cafe_large_mode per_family: no lambda given, so each family estimates its own.
    # A lambda estimated jointly over all large families didn't converge (CAFE5
    # discussion #132), which is why these families have always run one at a time.
    cafe5 \\
        -i ${hog_counts} \\
        -t ${species_tree} \\
        --cores ${task.cpus} \\
        ${e_flag} \\
        ${l_flag} \\
        ${z_flag} \\
        ${P_flag} \\
        -o ${outdir} \\
        2>&1 | tee ${logf} || true

    # What became of this family, for CAFE_FAMILY_STATUS. A family whose task fails
    # outright writes nothing and is reported there as task_failed.
    hog=\$(awk -F'\t' 'NR == 2 { print \$2 }' ${hog_counts})
    if grep -q "not present at the root from: 1 to 0" ${logf}; then
        echo "WARNING: ${hog_counts.baseName} has no genes on one side of the root, so CAFE5 dropped it (--cafe_zero_root keeps such families) — not in the large-family merge." >&2
        rm -rf ${outdir}
        status=not_present_at_root
    elif [ ! -f ${outdir}/Base_results.txt ] || \\
       ! grep -q "Final Likelihood" ${outdir}/Base_results.txt || \\
       grep -q "inf" ${outdir}/Base_results.txt; then
        echo "WARNING: ${hog_counts.baseName} has no finite likelihood${fixed_lambda ? " at the fixed lambda ${fixed_lambda}" : " even with its own lambda"} — dropping it from the large-family merge." >&2
        rm -rf ${outdir}
        status=did_not_converge
    else
        status=fitted
    fi
    printf "%s\\t%s\\n" "\$hog" "\$status" > large_status_${hog_counts.baseName}.tsv
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
    printf "%s\\tfitted\\n" "${hog_counts.baseName}" > large_status_${hog_counts.baseName}.tsv
    """
}
