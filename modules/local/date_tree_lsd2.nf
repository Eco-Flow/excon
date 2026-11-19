// --dating_method lsd2: time-calibrate the IQ-TREE species tree with IQ-TREE's
// least-square dating (LSD2) instead of ape::chronos, from the same calibration table.
// Three steps, because the IQ-TREE container has no R: check the calibrations and write
// IQ-TREE's date file; date; check the dated tree. See bin/lsd2_dating.R.

process LSD2_PREPARE {
    label 'process_single'
    container 'ecoflowucl/cafe:r-4.3.1'

    input:
    path tree_newick
    path calibrations

    output:
    path 'lsd2_dates.txt',    emit: dates
    path 'lsd2_outgroup.txt', emit: outgroup

    script:
    """
    Rscript ${projectDir}/bin/lsd2_dating.R prepare ${tree_newick} ${calibrations}
    """
}

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
}

process LSD2_CHECK {
    label 'process_single'
    container 'ecoflowucl/cafe:r-4.3.1'

    input:
    path nex
    path calibrations
    path tree_newick
    path iqtree_log

    output:
    path 'SpeciesTree_dated.nwk'   , emit: dated_tree
    path 'dating_calibrations.tsv' , emit: calibrations
    path 'dating_qc.tsv'           , emit: qc

    script:
    """
    Rscript ${projectDir}/bin/lsd2_dating.R check ${nex} ${calibrations} ${tree_newick} ${iqtree_log}
    """
}
