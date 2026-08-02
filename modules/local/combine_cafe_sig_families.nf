process COMBINE_CAFE_SIG_FAMILIES {
    label 'process_single'
    container 'ecoflowucl/cafe:r-4.3.1'

    input:
    path main_changes,  stageAs: 'main_changes.tsv'
    path main_sig,      stageAs: 'main_sig_changes.tsv'
    path large_changes, stageAs: 'large_changes.tsv'
    path large_sig,     stageAs: 'large_sig_changes.tsv'

    output:
    path "combined_changes_per_node.tsv",             emit: changes
    path "combined_significant_changes_per_node.tsv", emit: sig_changes

    script:
    """
    ${projectDir}/bin/combine_cafe_sig_families.py \\
        --main-changes main_changes.tsv \\
        --large-changes large_changes.tsv \\
        --main-sig-changes main_sig_changes.tsv \\
        --large-sig-changes large_sig_changes.tsv \\
        --out-changes combined_changes_per_node.tsv \\
        --out-sig-changes combined_significant_changes_per_node.tsv
    """

    stub:
    """
    touch combined_changes_per_node.tsv
    touch combined_significant_changes_per_node.tsv
    """
}
