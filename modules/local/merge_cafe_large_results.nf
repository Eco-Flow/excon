process MERGE_CAFE_LARGE_RESULTS {
    label 'process_single'
    container 'ecoflowucl/cafe:r-4.3.1'

    // The large-family track is best-effort supplementary analysis (as it was
    // before this rewrite) — an unexpected merge failure here should not sink
    // the main CAFE run or GO enrichment. Downstream _LARGE steps simply see
    // empty channels and are skipped, same as when nothing converged at all.
    errorStrategy 'ignore'

    input:
    path family_dirs   // collected Out_cafe_large_<HOG>/ dirs, one per converged family

    output:
    path "Out_cafe_large/",                 emit: results
    path "large_family_lambda_summary.tsv", emit: lambda_summary, optional: true
    path "converged.txt", optional: true,   emit: converged

    script:
    """
    ${projectDir}/bin/merge_cafe_large_results.py \\
        --outdir Out_cafe_large \\
        --lambda-summary large_family_lambda_summary.tsv \\
        Out_cafe_large_*/

    if [ -f Out_cafe_large/Base_asr.tre ]; then
        touch converged.txt
    fi
    """

    stub:
    """
    mkdir -p Out_cafe_large
    touch Out_cafe_large/Base_asr.tre
    touch Out_cafe_large/Base_branch_probabilities.tab
    touch Out_cafe_large/Base_change.tab
    touch Out_cafe_large/Base_count.tab
    touch Out_cafe_large/Base_family_results.txt
    touch Out_cafe_large/Base_clade_results.txt
    touch large_family_lambda_summary.tsv
    touch converged.txt
    """
}
