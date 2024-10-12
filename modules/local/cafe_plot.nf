process CAFE_PLOT {
    label 'process_single'
    container = 'ecoflowucl/cafeplotter:latest'

    input:
    path Cafe_dir

    output:
    path("cafe_plotter") , emit: cafe_plot_results
    path "versions.yml", emit: versions

    script:
    """
    cafeplotter -i ${Cafe_dir} -o cafe_plotter --format 'pdf'

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        Python version: \$(perl --version | grep "version" | sed 's/.*(//g' | sed 's/[)].*//')
    END_VERSIONS
    """
}
