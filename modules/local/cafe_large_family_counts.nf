process CAFE_LARGE_FAMILY_COUNTS {
    label 'process_single'
    container 'ecoflowucl/cafe:r-4.3.1'

    // A separate step rather than part of CAFE_PREP, so it can change without
    // invalidating CAFE_PREP's cache, and with it every CAFE5 run downstream.
    input:
    path n0_table           // N0.tsv as CAFE_PREP used it
    path filtering_report   // hog_filtering_report.tsv: why each family was left out
    path cafe_tree          // cafe_input_tree.txt: its tips are the count columns

    output:
    path "hog_gene_counts_large.tsv", emit: large_counts, optional: true
    path "large_families.txt",        emit: summary
    tuple val("${task.process}"), val('R'), val('4.3.1'), emit: versions_R, topic: versions

    script:
    """
    ${projectDir}/bin/cafe_large_families.R ${n0_table} ${filtering_report} ${cafe_tree}
    """
}
