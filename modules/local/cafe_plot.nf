process CAFE_PLOT {
    label 'process_low'
    container = 'ecoflowucl/cafeplotter:latest'

    input:
    path Cafe_dir

    output:
    path("*.pdf") , emit: cafe_go_pdfs
    path "versions.yml", emit: versions

    script:
    """
    cafeplotter -i ${Cafe_dir} -o cafe_plotter

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        Python version: \$(perl --version | grep "version" | sed 's/.*(//g' | sed 's/[)].*//')
    END_VERSIONS
    """
}
