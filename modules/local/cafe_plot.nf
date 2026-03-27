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
    cafeplotter -i ${Cafe_dir} -o cafe_plotter --format 'pdf'
    """
}
