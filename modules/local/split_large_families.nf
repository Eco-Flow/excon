process SPLIT_LARGE_FAMILIES {
    label 'process_single'
    container 'ecoflowucl/cafe:r-4.3.1'

    input:
    path large_counts   // hog_gene_counts_large.tsv from CAFE_PREP (optional — may be absent)

    output:
    path "large_family_splits/*.tsv", emit: family_tables, optional: true

    script:
    """
    ${projectDir}/bin/split_large_families.py -o large_family_splits ${large_counts}
    """

    stub:
    """
    mkdir -p large_family_splits
    printf "Desc\\tHOG\\tspA\\tspB\\n" > large_family_splits/N0.HOG0000001.tsv
    printf "n/a\\tN0.HOG0000001\\t0\\t99\\n" >> large_family_splits/N0.HOG0000001.tsv
    """
}
