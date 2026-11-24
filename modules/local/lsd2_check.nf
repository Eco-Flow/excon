// --dating_method lsd2, step 3 of 3 (subworkflows/local/lsd2_dating): check the
// dated tree (ultrametric, same topology and root, calibrations held) and write
// the same outputs as DATE_TREE (ape::chronos).

process LSD2_CHECK {
    label 'process_single'

    container 'ecoflowucl/cafe:r-4.3.1' //R version 4.3.1, includes ape

    input:
    path nex
    path calibrations
    path tree_newick
    path iqtree_log

    output:
    path 'SpeciesTree_dated.nwk'   , emit: dated_tree
    path 'dating_calibrations.tsv' , emit: calibrations
    path 'dating_qc.tsv'           , emit: qc
    tuple val("${task.process}"), val('ape'), eval("Rscript -e 'cat(as.character(packageVersion(\"ape\")))'"), emit: versions, topic: versions

    when:
    task.ext.when == null || task.ext.when

    script:
    """
    Rscript ${projectDir}/bin/lsd2_dating.R check ${nex} ${calibrations} ${tree_newick} ${iqtree_log}
    """

    stub:
    """
    touch SpeciesTree_dated.nwk
    touch dating_calibrations.tsv
    touch dating_qc.tsv
    """
}
