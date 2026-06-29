process SUMMARIZE_CAFE_GO {
    tag "$tag"
    label 'process_single'
    container 'ecoflowucl/chopgo:r-4.3.2_python-3.10_perl-5.38'

    input:
    tuple val(tag), path(topgo_files)

    output:
    tuple val(tag), path("Go_summary_pos.tsv"),          emit: pos_tsv
    tuple val(tag), path("Go_summary_neg.tsv"),          emit: neg_tsv
    tuple val(tag), path("Go_summary_posneg_merged.tsv"), emit: merged_tsv
    tuple val("${task.process}"), val('perl'), eval("perl --version 2>&1 | grep 'version' | sed 's/.*(//; s/[)].*//'"), emit: versions_perl, topic: versions

    script:
    """
    sum_cafe.pl
    """
}
