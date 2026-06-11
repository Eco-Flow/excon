process CAFE_NODE_GUIDE {
    label 'process_single'
    container 'ecoflowucl/cafe:r-4.3.1'

    input:
    path cafe_dir

    output:
    path("cafe_node_label_guide.pdf"), emit: node_guide, optional: true
    tuple val("${task.process}"), val('R'), eval("R --version 2>&1 | grep 'R version' | sed 's/R version \\([0-9.]*\\).*/\\1/'"), emit: versions_R, topic: versions

    script:
    """
    ${projectDir}/bin/cafe_node_label_guide.R ${cafe_dir}
    """
}
