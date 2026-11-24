// --dating_method lsd2, step 1 of 3 (subworkflows/local/lsd2_dating): check the
// calibrations against the rooted tree and write IQ-TREE's --date file and the
// outgroup for -o. Separate from DATE_TREE_LSD2 because the IQ-TREE container has no R.

process LSD2_PREPARE {
    label 'process_single'

    container 'ecoflowucl/cafe:r-4.3.1' //R version 4.3.1, includes ape

    input:
    path tree_newick
    path calibrations

    output:
    path 'lsd2_dates.txt'    , emit: dates
    path 'lsd2_outgroup.txt' , emit: outgroup
    tuple val("${task.process}"), val('ape'), eval("Rscript -e 'cat(as.character(packageVersion(\"ape\")))'"), emit: versions, topic: versions

    when:
    task.ext.when == null || task.ext.when

    script:
    """
    Rscript ${projectDir}/bin/lsd2_dating.R prepare ${tree_newick} ${calibrations}
    """

    stub:
    """
    touch lsd2_dates.txt
    touch lsd2_outgroup.txt
    """
}
