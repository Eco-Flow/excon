process CAFE_MODEL_COMPARE {
    tag "cafe_model_compare"
    label 'process_single'
    container 'ecoflowucl/cafe:r-4.3.1'

    input:
    path uniform_results   // Out_cafe_kN       from CAFE_RUN_K (best k, no Poisson)
    path poisson_results   // Out_cafe_kN_poisson from CAFE_RUN_BEST

    output:
    path "cafe_model_comparison.tsv", emit: comparison_table
    path "best_model.txt",            emit: best_model
    path "best_cafe_results/",        emit: best_results
    tuple val("${task.process}"), val('cafe'), val('4.2.1'), emit: versions_cafe, topic: versions

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

    # Choose the best model. Lower -lnL = better fit. Pick whichever has a usable
    # score; only fall over if neither does.
    if is_finite "\$uniform_score" && is_finite "\$poisson_score"; then
        best=\$(awk -v u="\$uniform_score" -v p="\$poisson_score" \\
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

    # NA only for table display when a score was missing/inf
    uniform_score=\${uniform_score:-NA}
    poisson_score=\${poisson_score:-NA}

    # Write comparison table
    printf "model\tdirectory\tneg_lnL\tSelected\n"          > cafe_model_comparison.tsv
    printf "uniform\t${uniform_results}\t\$uniform_score\t"  >> cafe_model_comparison.tsv
    [ "\$best" = "uniform" ] && printf "BEST\n" >> cafe_model_comparison.tsv || printf "\n" >> cafe_model_comparison.tsv
    printf "poisson\t${poisson_results}\t\$poisson_score\t"  >> cafe_model_comparison.tsv
    [ "\$best" = "poisson" ] && printf "BEST\n" >> cafe_model_comparison.tsv || printf "\n" >> cafe_model_comparison.tsv

    echo "Model comparison complete — best model: \$best (-lnL: uniform=\$uniform_score, poisson=\$poisson_score)" >&2

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
