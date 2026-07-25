process CAFE_SIG_FAMILIES {
    label 'process_single'
    container 'ecoflowucl/cafe:r-4.3.1'

    input:
    path cafe_results
    path n0_table
    // Both are optional and default to the same assets/NO_FILE placeholder, so they
    // must be staged under distinct names or Nextflow rejects the name collision.
    // cafe_sig_families.R guards with dir.exists(), and the placeholder is a file,
    // so an unset input is still correctly treated as absent.
    path msa_dir,  stageAs: 'og_msa_dir'
    path tree_dir, stageAs: 'og_gene_tree_dir'

    output:
    path("changes_per_node.tsv"),             emit: changes
    path("significant_changes_per_node.tsv"), emit: sig_changes
    path("focus_families_manifest.tsv"),      emit: focus_manifest, optional: true
    path("focus_families/**"),                emit: focus_families,  optional: true
    tuple val("${task.process}"), val('R'), val('4.3.1'), emit: versions_R, topic: versions

    script:
    def focus = params.cafe_focus_clades ? "\"${params.cafe_focus_clades}\"" : "NA"
    def pcut  = params.go_cutoff ?: 0.05
    """
    ${projectDir}/bin/cafe_sig_families.R \\
        ${cafe_results} \\
        ${n0_table} \\
        ${pcut} \\
        ${focus} \\
        ${msa_dir} \\
        ${tree_dir}
    """
}
