process CAFE_PREP {
    label 'process_high'
    label 'process_long'
    container 'ecoflowucl/cafe:r-4.3.1' //R version 4.3.1, cafe version 4.2.1 (confusing)


    errorStrategy {
        if (task.attempt <= 3 && task.exitStatus == 1) {
            log.warn "CAFE_PREP: convergence failed on attempt ${task.attempt} — retrying with stricter differential filtering"
            return 'retry'
        }
        return 'ignore'
    }
    maxRetries 3

    input:
    path table
    path tree_newick

    output:
    path("hog_gene_counts.tsv"),                         emit: prepared_counts
    path("cafe_input_tree.txt"),                         emit: cafe_tree
    path("SpeciesTree_rooted_ultra.txt"),                emit: prepared_tree
    path("pruned_tree"),                                 emit: pruned_tree
    path("N0.tsv"),                                      emit: N0_table
    path("Out_cafe"),                                    emit: results
    path("Out_cafe/Base_count.tab"),                     emit: result_nftest
    path("Out_cafe_errormodel/Base_error_model.txt"),    emit: error_model
    path("hog_filtering_report.tsv"),                    emit: filtering_report
    path("hog_gene_counts_large.tsv"),                   emit: large_counts,     optional: true
    path("lambda.txt"),                                  emit: lambda
    path("cafe_base.log"),                               emit: base_log
    path("cafe_errormodel.log"),                         emit: errormodel_log
    tuple val("${task.process}"), val('R'),    val('4.3.1'), emit: versions_R,    topic: versions
    tuple val("${task.process}"), val('cafe'), val('4.2.1'), emit: versions_cafe, topic: versions

    script:
    def base_differential = params.cafe_max_differential ?: 50
    // By default the first attempt is unfiltered, so no family is discarded when
    // CAFE5 can cope with the full set; each retry then halves the threshold.
    // --cafe_filter_first starts filtering immediately, which is worth setting when
    // the data are already known to need it — the unfiltered attempt is otherwise a
    // guaranteed failure, and relying on retries is fragile if a run gets interrupted.
    def first_attempt_no = params.cafe_filter_first ? 1 : 2
    def max_differential  = (base_differential / Math.pow(2, task.attempt - first_attempt_no)).toInteger()
    def use_filtering    = params.cafe_filter_first || task.attempt > 1
    // A tree from DATE_TREE is already ultrametric and in Myr, so it is handled
    // exactly like a user-supplied dated tree.
    def is_dated         = (params.input_tree_is_dated || params.tree_calibrations) ? 'true' : 'false'
    def z_flag           = params.cafe_zero_root ? '-z' : ''
    // Matches maxRetries below: the final attempt must not exit 1, or the run ends
    // with no result at all rather than one lacking an error model.
    def can_retry        = task.attempt <= 3
    """
    export PATH=\$PATH:/usr/bin
    set -e

    [ "${table}" = "N0.tsv" ] || cp ${table} N0.tsv
    cp ${tree_newick} pruned_tree
    sed -i 's/\\.clean//g' pruned_tree
    sed -i 's/\\.clean//g' N0.tsv

    # Scale factor 1: RESCALE_TREE has already applied --tree_scale_factor to the
    # incoming tree, and chronoMPL() is linear, so applying it again here would
    # scale the tree by the factor squared.
    if [ "${use_filtering}" = "true" ]; then
        echo "CAFE_PREP attempt ${task.attempt}: applying differential filtering (threshold: ${max_differential})"
        Rscript ${projectDir}/bin/cafe_prep_filtered.R ${max_differential} 1 ${is_dated}
    else
        echo "CAFE_PREP attempt ${task.attempt}: no filtering"
        Rscript ${projectDir}/bin/cafe_prep.R 1 ${is_dated}
    fi

    # ---------------------------------------------------------------
    # Stage 1: base run (λ estimation, no error model)
    # Uses cafe_input_tree.txt — the ultrametric tree cafe_prep.R emits.
    # When --input_tree_is_dated, this is the supplied time-calibrated tree
    # UNCHANGED (branch lengths in Myr); otherwise it is the chronoMPL
    # ultrametric tree scaled once by --tree_scale_factor. Its tips match
    # the (subset) gene-count columns exactly.
    # ---------------------------------------------------------------
    cafe5 \\
        -i hog_gene_counts.tsv \\
        -t cafe_input_tree.txt \\
        --cores ${task.cpus} \\
        ${z_flag} \\
        -o Out_cafe \\
        2>&1 | tee cafe_base.log
    cafe5_exit=\${PIPESTATUS[0]}
    [ \$cafe5_exit -ne 0 ] && exit \$cafe5_exit


    # ---------------------------------------------------------------
    # Check for convergence failure — CAFE5 exits 0 even on -inf lnL
    # These conditions trigger a retry with differential filtering
    # ---------------------------------------------------------------
    if grep -q "largest size differential" cafe_base.log; then
        echo "ERROR: CAFE5 detected size differential error — retrying with filtering" >&2
        exit 1
    fi

    # Check for actual -inf likelihood values, not the word "Inferring"
    if grep -qE "Score \\(-lnL\\):\\s+inf" cafe_base.log && ! grep -q "Final -lnL:" cafe_base.log; then
        echo "ERROR: CAFE5 failed to converge (infinite likelihoods) — retrying with filtering" >&2
        exit 1
    fi

    if ! grep -q "Final -lnL:" cafe_base.log; then
        echo "ERROR: CAFE5 produced no likelihood score — retrying with filtering" >&2
        exit 1
    fi

    # Extract lambda estimate for fixed-lambda re-analysis of large families
    grep "^Lambda:" Out_cafe/Base_results.txt | awk '{print \$2}' > lambda.txt


    # ---------------------------------------------------------------
    # Stage 2: estimate error model
    # Quantifies assembly/annotation error in gene family counts.
    # The resulting Base_error_model.txt is passed to CAFE_RUN_K so
    # that all downstream k-sweep runs correct for this error.
    # ---------------------------------------------------------------
    # set +e so a failing error model does not abort the script: process.shell sets
    # -e and pipefail, which would otherwise kill the task at the pipeline below and
    # skip the fallback entirely.
    set +e
    cafe5 \\
        -i hog_gene_counts.tsv \\
        -t cafe_input_tree.txt \\
        --cores ${task.cpus} \\
        ${z_flag} \\
        -e \\
        -o Out_cafe_errormodel \\
        2>&1 | tee cafe_errormodel.log
    errormodel_exit=\${PIPESTATUS[0]}
    set -e

    # A failed error model is non-fatal — downstream processes handle an empty file
    # via the optional NO_FILE fallback pattern. CAFE5 can also exit 0 without writing
    # the model at all, so the file itself is checked rather than just the exit status.
    touch cafe_errormodel.log
    if [ \$errormodel_exit -ne 0 ] || [ ! -s Out_cafe_errormodel/Base_error_model.txt ]; then
        if [ "${can_retry}" = "true" ]; then
            # CAFE5 exits 0 even when it fails to converge here, so the retry has to be
            # triggered explicitly, exactly as the base run does above. Without this the
            # task exits 0, Nextflow fails it for a missing output, and errorStrategy
            # sees exitStatus 0 rather than 1 and ignores it instead of retrying.
            echo "ERROR: error model did not converge — retrying with stricter differential filtering" >&2
            exit 1
        fi
        echo "WARNING: no error model produced after ${task.attempt} attempts — continuing without it" >&2
        mkdir -p Out_cafe_errormodel
        touch Out_cafe_errormodel/Base_error_model.txt
    fi

    exit 0
    """
}
