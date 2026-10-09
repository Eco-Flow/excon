process CAFE_MODEL_COMPARE {
    tag "cafe_model_compare"
    label 'process_single'
    container 'ecoflowucl/cafe:r-4.3.1'

    input:
    path uniform_results   // Out_cafe_kN       from CAFE_RUN_K (best k, no Poisson)
    path poisson_results   // Out_cafe_kN_poisson from CAFE_RUN_BEST
    path prepared_counts   // hog_gene_counts.tsv from CAFE_PREP: the families CAFE5 was given

    output:
    path "cafe_model_comparison.tsv", emit: comparison_table
    path "best_model.txt",            emit: best_model
    path "best_cafe_results/",        emit: best_results
    path "families_not_modelled.tsv", emit: not_modelled, optional: true
    path "*_warning.txt",             emit: warnings,     optional: true
    tuple val("${task.process}"), val('cafe5'), val('1.1.0'), emit: versions_cafe, topic: versions

    script:
    """
    # Parse -lnL from a CAFE5 results directory.
    # Line looks like: "Model Gamma Final Likelihood (-lnL): 185019.11912528"
    # The value may be an integer, a decimal, or "inf" (failed to converge), so
    # take everything after "-lnL):" rather than assuming a decimal point.
    parse_score() {
        grep -h "Final Likelihood" "\$1"/Base_results.txt "\$1"/Gamma_results.txt 2>/dev/null \\
            | head -1 \\
            | sed -E 's/.*-lnL\\):[[:space:]]*//' \\
            | grep -oE '^(inf|[0-9]+(\\.[0-9]+)?)' \\
            || true
    }

    # A score is usable only if it is a finite number (not empty, not "inf"/"nan").
    is_finite() { case "\$1" in ''|inf|-inf|nan) return 1 ;; *) return 0 ;; esac ; }

    uniform_score=\$(parse_score ${uniform_results})
    poisson_score=\$(parse_score ${poisson_results})

    # Parameters: lambda, plus alpha when k > 1 (gamma), plus the root-size mean for
    # the Poisson root. The error model is read from a file, not estimated, so it adds
    # none. AIC/BIC put the two fits on the same footing; -lnL alone favours the extra
    # parameter, and with thousands of families the root-size prior can make up most of
    # the difference between them (see root_split_warning.txt when it is written).
    k=\$(echo "${uniform_results}" | sed -E 's/.*_k([0-9]+).*/\\1/')
    n_uniform=\$([ "\$k" -gt 1 ] && echo 2 || echo 1)
    n_poisson=\$((n_uniform + 1))
    n_families=\$(cat ${uniform_results}/*_family_results.txt 2>/dev/null | grep -vc '^#' || true)

    # AIC = 2n + 2(-lnL); BIC = n ln(families) + 2(-lnL); "NA" when the fit failed
    aic() { is_finite "\$1" && awk -v s="\$1" -v n="\$2" 'BEGIN { printf "%.4f", 2*n + 2*s }' || echo NA ; }
    bic() { is_finite "\$1" && [ "\${n_families:-0}" -gt 0 ] && awk -v s="\$1" -v n="\$2" -v f="\$n_families" 'BEGIN { printf "%.4f", n*log(f) + 2*s }' || echo NA ; }
    uniform_aic=\$(aic "\$uniform_score" \$n_uniform)
    poisson_aic=\$(aic "\$poisson_score" \$n_poisson)

    # Choose the best model by AIC, as for k. Pick whichever has a usable score;
    # only fall over if neither does, so a failed fit (NA) is never selected.
    if is_finite "\$uniform_score" && is_finite "\$poisson_score"; then
        best=\$(awk -v u="\$uniform_aic" -v p="\$poisson_aic" \\
            'BEGIN { print (p+0 < u+0) ? "poisson" : "uniform" }')
    elif is_finite "\$uniform_score"; then
        best=uniform
    elif is_finite "\$poisson_score"; then
        best=poisson
    else
        echo "ERROR: neither model produced a usable -lnL (uniform='\$uniform_score' poisson='\$poisson_score')." >&2
        echo "       The CAFE k-runs feeding CAFE_MODEL_COMPARE are empty or failed to converge." >&2
        exit 1
    fi
    echo "\$best" > best_model.txt

    # Write comparison table
    printf "model\tdirectory\tk\tn_params\tneg_lnL\tAIC\tBIC\tn_families\tSelected\n" > cafe_model_comparison.tsv
    printf "uniform\t${uniform_results}\t\$k\t\$n_uniform\t\${uniform_score:-NA}\t\$uniform_aic\t\$(bic "\$uniform_score" \$n_uniform)\t\${n_families:-NA}\t\$([ "\$best" = uniform ] && echo BEST)\n" >> cafe_model_comparison.tsv
    printf "poisson\t${poisson_results}\t\$k\t\$n_poisson\t\${poisson_score:-NA}\t\$poisson_aic\t\$(bic "\$poisson_score" \$n_poisson)\t\${n_families:-NA}\t\$([ "\$best" = poisson ] && echo BEST)\n" >> cafe_model_comparison.tsv

    echo "Model comparison complete — best model: \$best (AIC: uniform=\$uniform_aic, poisson=\$poisson_aic)" >&2

    # Copy the winning model's results directory for publishing
    if [ "\$best" = "uniform" ]; then
        src=${uniform_results}
    else
        src=${poisson_results}
    fi
    cp -rL "\$src" best_cafe_results

    # Guard: the published best/ must contain a real CAFE result. An empty dir
    # here means the staged input was empty (e.g. a stale -resume cache pointing
    # at a CAFE run that hadn't populated yet) — fail loudly instead of silently
    # publishing an empty best/.
    if ! grep -q "Final Likelihood" best_cafe_results/Base_results.txt best_cafe_results/Gamma_results.txt 2>/dev/null; then
        echo "ERROR: chosen model '\$best' (\$src) produced no usable results;" >&2
        echo "       best_cafe_results has no *_results.txt with a likelihood." >&2
        echo "       This usually means a stale -resume cache staged an empty CAFE run dir;" >&2
        echo "       remove the cached CAFE_MODEL_COMPARE work dir and re-run." >&2
        ls -la best_cafe_results >&2 || true
        exit 1
    fi

    # Writes root_split_warning.txt when one species is alone on one side of the root
    ${projectDir}/bin/cafe_root_check.R best_cafe_results

    # Families CAFE5 was given but left out of its results: those not present at the
    # root (no genes on one side of it), which it drops unless run with -z. Its log says
    # how many; this says which, as hog_filtering_report.tsv still counts them as retained.
    tail -n +2 ${prepared_counts} | cut -f2 | sort > input_families.txt
    cat best_cafe_results/*_change.tab | tail -n +2 | cut -f1 | sort > modelled_families.txt
    comm -23 input_families.txt modelled_families.txt > not_modelled.txt
    n_dropped=\$(wc -l < not_modelled.txt | tr -d ' ')
    if [ "\$n_dropped" -gt 0 ]; then
        { printf "HOG\\treason\\n"; awk '{ print \$1 "\\tnot present at the root (no genes on one side of it)" }' not_modelled.txt; } > families_not_modelled.tsv
        printf "WARNING: CAFE5 left out %s of the %s families it was given because they have no genes on one side\\nof the root, so it can't place them there. They are listed in cafe/model_comparison/families_not_modelled.tsv.\\nRun with --cafe_zero_root to keep them (CAFE5's -z).\\n" \\
            "\$n_dropped" "\$(wc -l < input_families.txt | tr -d ' ')" > root_filter_warning.txt
    fi
    rm -f input_families.txt modelled_families.txt not_modelled.txt
    """

    stub:
    """
    echo -e "model\tdirectory\tneg_lnL\tSelected" > cafe_model_comparison.tsv
    echo -e "uniform\tOut_cafe_k3\t195812.3\tBEST"  >> cafe_model_comparison.tsv
    echo -e "poisson\tOut_cafe_k3_poisson\t196100.1\t" >> cafe_model_comparison.tsv
    echo "uniform" > best_model.txt
    mkdir -p best_cafe_results
    touch best_cafe_results/Base_results.txt
    """
}
