process SPLIT_LARGE_FAMILIES {
    label 'process_single'

    input:
    path large_counts   // hog_gene_counts_large.tsv from CAFE_PREP (optional — may be absent)

    output:
    path "large_family_splits/*.tsv", emit: family_tables, optional: true

    script:
    """
    mkdir -p large_family_splits
    awk -F'\\t' '
        NR == 1 { header = \$0; next }
        {
            hog = \$2
            gsub(/[^A-Za-z0-9_.-]/, "_", hog)
            fname = "large_family_splits/" hog ".tsv"
            print header > fname
            print \$0 >> fname
            close(fname)
        }
    ' ${large_counts}
    """

    stub:
    """
    mkdir -p large_family_splits
    printf "Desc\\tHOG\\tspA\\tspB\\n" > large_family_splits/N0.HOG0000001.tsv
    printf "n/a\\tN0.HOG0000001\\t0\\t99\\n" >> large_family_splits/N0.HOG0000001.tsv
    """
}
