process CAFE_GO_PREP {
    tag "cafe_go_prep"
    label 'process_low'
    container 'ecoflowucl/chopgo:r-4.3.2_python-3.10_perl-5.38'

    input:
    path cafe_results
    path n0_table
    path go

    output:
    path("chopgo_jobs.tsv"),    emit: manifest
    path("*.pos.txt"),          emit: pos_files
    path("*.neg.txt"),          emit: neg_files
    path("*.BK.txt.uniq"),      emit: bk_files
    path("OG_GO_format.tsv"),   emit: og_go
    path("CAFE_summary.txt"),   emit: cafe_summary
    tuple val("${task.process}"), val('perl'), eval("perl --version 2>&1 | grep 'version' | sed 's/.*(//; s/[)].*//'"), emit: versions_cafeprep_perl, topic: versions

    script:
    """
    ln -s ${cafe_results} Out_cafe
    ${projectDir}/bin/cafe_go_prep.pl ${params.go_cutoff} ${params.go_type} ${params.go_max_plot}
    """
}
