// --dating_method lsd2, step 2 of 3 (subworkflows/local/lsd2_dating): time-calibrate
// the IQ-TREE species tree with IQ-TREE's least-square dating (LSD2), from the
// supermatrix it was inferred from and the dates LSD2_PREPARE wrote.

process DATE_TREE_LSD2 {
    label 'process_low'

    container "${ workflow.containerEngine in ['singularity', 'apptainer'] && !task.ext.singularity_pull_docker_container ?
        'https://depot.galaxyproject.org/singularity/iqtree:3.1.3--h8471819_0' :
        'quay.io/biocontainers/iqtree:3.1.3--h8471819_0' }"

    input:
    path alignment      // the supermatrix the tree was inferred from
    path tree_newick    // the rooted species tree (substitutions/site)
    path dates
    path outgroup

    output:
    path 'dated.timetree.nex', emit: nex
    path 'dated.timetree.lsd', emit: report
    path 'dated.log',          emit: log
    tuple val("${task.process}"), val('iqtree'), eval("iqtree -version 2>&1 | grep -oE '[0-9]+[.][0-9]+[.][0-9]+' | head -n1"), emit: versions, topic: versions

    when:
    task.ext.when == null || task.ext.when

    script:
    // -blfix keeps the tree's branch lengths, so the model only serves LSD2's variances;
    // -o names one side of the root, so LSD2 keeps the root rather than estimating it
    """
    iqtree \\
        -s ${alignment} \\
        -te ${tree_newick} \\
        -blfix \\
        -m LG \\
        --date ${dates} \\
        --date-tip 0 \\
        -o "\$(cat ${outgroup})" \\
        -T ${task.cpus} \\
        -pre dated
    """

    stub:
    """
    touch dated.timetree.nex
    touch dated.timetree.lsd
    touch dated.log
    """
}
