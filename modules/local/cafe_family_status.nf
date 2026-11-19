process CAFE_FAMILY_STATUS {
    tag "cafe_family_status"
    label 'process_single'
    container 'ecoflowucl/chopgo:r-4.3.2_python-3.10_perl-5.38'

    input:
    path filtering_report                       // hog_filtering_report.tsv from CAFE_PREP
    path best_results                           // the selected main model's CAFE5 results
    path large_counts                           // hog_gene_counts_large.tsv, or NO_FILE
    path large_status, stageAs: 'large_status/*' // large_status_*.tsv from CAFE_RUN_LARGE (may be none)

    output:
    path "cafe_family_status.tsv",         emit: status
    path "cafe_family_status_summary.txt", emit: summary
    tuple val("${task.process}"), val('python'), val('3.10.13'), emit: versions_python, topic: versions

    script:
    def large  = large_counts.name != 'NO_FILE' ? "--large-counts ${large_counts}" : ''
    """
    /usr/bin/python3.10 ${projectDir}/bin/cafe_family_status.py \\
        --filtering-report ${filtering_report} \\
        --best-results ${best_results} \\
        ${large} \\
        --large-status \$(find -L large_status -name '*.tsv' 2>/dev/null) \\
        -o cafe_family_status.tsv
    """

    stub:
    """
    printf "HOG\\tstatus\\ttrack\\tprep_filter\\tn_species_present\\tn_max\\tn_min\\ttotal_genes\\tdifferential\\n" > cafe_family_status.tsv
    echo "CAFE family status: 0 families in the orthogroup table" > cafe_family_status_summary.txt
    """
}
