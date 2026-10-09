process CAFE_PLOT {
    label 'process_single'
    container 'ecoflowucl/cafeplotter:latest'

    input:
    path Cafe_dir

    output:
    path("cafe_plotter") , emit: cafe_plot_results
    tuple val("${task.process}"), val('cafeplotter'), eval("cafeplotter --version"), emit: versions_cafeplotter, topic: versions

    script:
    """
    # CAFE5 estimates each family's p-value from 1,000 simulated families, so a p-value
    # of 0 means p < 0.001; cafeplotter printed it in the title as "p-value=0". Patch
    # its display through a sitecustomize.py on PYTHONPATH (its install is read-only).
    mkdir pypatch
    cat > pypatch/sitecustomize.py <<'PYEOF'
from cafeplotter.cafeparser import FamilyResult
from cafeplotter.treeplotter import TreePlotter
_display = FamilyResult.display_pvalue.fget
FamilyResult.display_pvalue = property(lambda self: "< 0.001" if self.pvalue == 0 else _display(self))
_set_title = TreePlotter.set_title
TreePlotter.set_title = lambda self, label, **kw: _set_title(self, label.replace("p-value=< ", "p-value < "), **kw)
PYEOF
    export PYTHONPATH="\$PWD/pypatch\${PYTHONPATH:+:\$PYTHONPATH}"
    # Fail fast if the patch didn't load: Python only warns when sitecustomize.py fails
    python3 -c "from cafeplotter.cafeparser import FamilyResult; assert FamilyResult('x', 0, True).display_pvalue == '< 0.001'"

    if compgen -G "${Cafe_dir}/*_asr.tre" > /dev/null 2>&1; then
        cafeplotter -i ${Cafe_dir} -o cafe_plotter --format 'pdf'
        cafeplotter -i ${Cafe_dir} -o cafe_plotter --format 'svg'

        # cafeplotter writes one pdf and one svg per significant family — hundreds
        # of small files on a large analysis. Archived into one file; the summary
        # plot and result_summary.tsv (the ones actually worth browsing directly)
        # are left as they are.
        if [ -d cafe_plotter/gene_family ]; then
            tar -czf cafe_plotter/gene_family.tar.gz -C cafe_plotter gene_family
            rm -rf cafe_plotter/gene_family
        fi
    else
        echo "No *_asr.tre found in ${Cafe_dir} — CAFE5 did not converge, skipping plot."
        mkdir -p cafe_plotter
        touch cafe_plotter/no_convergence.txt
    fi
    """
}
