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
    if compgen -G "${Cafe_dir}/*_asr.tre" > /dev/null 2>&1; then
        cafeplotter -i ${Cafe_dir} -o cafe_plotter --format 'pdf'
    else
        echo "No *_asr.tre found in ${Cafe_dir} — CAFE5 did not converge, skipping plot."
        mkdir -p cafe_plotter
        touch cafe_plotter/no_convergence.txt
    fi
    """
}
