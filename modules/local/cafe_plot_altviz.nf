process CAFE_PLOT_ALTVIZ {
    label 'process_single'
    container 'ecoflowucl/cafe:r-4.3.1'

    input:
    path cafe_dir
    path cafe_summary

    output:
    path("cafe_sig_hog_tree.pdf"),  emit: hog_pdf,  optional: true
    path("cafe_sig_hog_tree.svg"),  emit: hog_svg,  optional: true
    path("cafe_sig_gene_tree.pdf"), emit: gene_pdf, optional: true
    path("cafe_sig_gene_tree.svg"), emit: gene_svg, optional: true
    tuple val("${task.process}"), val('R'), eval("R --version 2>&1 | grep 'R version' | sed 's/R version \\([0-9.]*\\).*/\\1/'"), emit: versions_R, topic: versions

    script:
    """
    ${projectDir}/bin/cafe_plot_altviz.R ${cafe_dir} ${cafe_summary}
    """
}
