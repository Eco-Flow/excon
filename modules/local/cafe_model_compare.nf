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
    tuple val("${task.process}"), val('cafe'), val('4.2.1'), emit: versions_cafe, topic: versions

    script:
    """
    # Parse -lnL from a CAFE5 results directory
    # Looks for: "Model * Final Likelihood (-lnL): X.X"
    parse_score() {
        grep -h "Final Likelihood" "\$1"/Base_results.txt "\$1"/Gamma_results.txt 2>/dev/null \\
            | head -1 \\
            | grep -oE '[0-9]+\\.[0-9]+' \\
            || true
    }

    uniform_score=\$(parse_score ${uniform_results})
    poisson_score=\$(parse_score ${poisson_results})

    if [ -z "\$uniform_score" ] || [ -z "\$poisson_score" ]; then
        echo "WARNING: Could not parse scores — defaulting to uniform model" >&2
        echo "uniform" > best_model.txt
        uniform_score=\${uniform_score:-NA}
        poisson_score=\${poisson_score:-NA}
    else
        # Lower -lnL = better fit
        best=\$(awk -v u="\$uniform_score" -v p="\$poisson_score" \\
            'BEGIN { print (p+0 < u+0) ? "poisson" : "uniform" }')
        echo "\$best" > best_model.txt
    fi

    best=\$(cat best_model.txt)

    # Write comparison table
    printf "model\tdirectory\tneg_lnL\tSelected\n"          > cafe_model_comparison.tsv
    printf "uniform\t${uniform_results}\t\$uniform_score\t"  >> cafe_model_comparison.tsv
    [ "\$best" = "uniform" ] && printf "BEST\n" >> cafe_model_comparison.tsv || printf "\n" >> cafe_model_comparison.tsv
    printf "poisson\t${poisson_results}\t\$poisson_score\t"  >> cafe_model_comparison.tsv
    [ "\$best" = "poisson" ] && printf "BEST\n" >> cafe_model_comparison.tsv || printf "\n" >> cafe_model_comparison.tsv

    echo "Model comparison complete — best model: \$best (-lnL: uniform=\$uniform_score, poisson=\$poisson_score)" >&2
    """

    stub:
    """
    echo -e "model\tdirectory\tneg_lnL\tSelected" > cafe_model_comparison.tsv
    echo -e "uniform\tOut_cafe_k3\t195812.3\tBEST"  >> cafe_model_comparison.tsv
    echo -e "poisson\tOut_cafe_k3_poisson\t196100.1\t" >> cafe_model_comparison.tsv
    echo "uniform" > best_model.txt
    """
}
