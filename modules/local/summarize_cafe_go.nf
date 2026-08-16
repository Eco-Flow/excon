process SUMMARIZE_CAFE_GO {
    tag "$tag"
    label 'process_single'
    container 'ecoflowucl/chopgo:r-4.3.2_python-3.10_perl-5.38'

    input:
    tuple val(tag), path(topgo_files)
    path pdf_files
    path svg_files

    output:
    tuple val(tag), path("Go_summary_pos.tsv"),          emit: pos_tsv
    tuple val(tag), path("Go_summary_neg.tsv"),          emit: neg_tsv
    tuple val(tag), path("Go_summary_posneg_merged.tsv"), emit: merged_tsv
    tuple val(tag), path("cafe_go_pos.tar.gz"), emit: pos_archive, optional: true
    tuple val(tag), path("cafe_go_neg.tar.gz"), emit: neg_archive, optional: true
    tuple val("${task.process}"), val('perl'), eval("perl --version 2>&1 | grep 'version' | sed 's/.*(//; s/[)].*//'"), emit: versions_perl, topic: versions

    script:
    """
    sum_cafe.pl ${params.go_algo}

    # One raw TopGO table plus a couple of barplots per species/direction adds up
    # to hundreds of loose files — every one of them carries '.pos.' or '.neg.'
    # in its name, so they sort cleanly into two archives without needing to
    # know their extensions up front.
    mkdir -p pos neg
    for f in *.pos.*; do [ -e "\$f" ] && cp "\$f" pos/; done
    for f in *.neg.*; do [ -e "\$f" ] && cp "\$f" neg/; done

    [ -n "\$(ls -A pos)" ] && tar -czf cafe_go_pos.tar.gz pos || true
    [ -n "\$(ls -A neg)" ] && tar -czf cafe_go_neg.tar.gz neg || true
    """
}
