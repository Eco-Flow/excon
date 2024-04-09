process GO_ASSIGN {
    label 'process_low'
    label 'error_ignore'
    tag "GO_assign_${Focal}"
    container = 'ecoflowucl/chopgo:r-4.3.2_python-3.10_perl-5.38'

    input:
    path Go_files
    path Orthogroups
    path Focal
    path Gene_to_trans

    output:
    path("${Focal}_Result_All_Combine_GO_format") , emit: go_hash
    path("OG_GO_format.tsv") , emit: go_og
    path("*results_ALL.tab.pdf") , emit: duplicate_go
    path("*family_expansions.txt") , emit: go_counts
    path("*transcripts_Combine_GO_format.txt"), emit: trans_go

    script:
    """
    perl -pe 's/\r\n|\n|\r/\n/g' ${Orthogroups} > Orthogroups.nomac.tsv
    ${projectDir}/bin/Goatee_ortho_go_match.pl Orthogroups.nomac.tsv ${Focal}
    ${projectDir}/bin/Orthofinder_duplicate_go.pl
    ${projectDir}/bin/Orthofinder_gene_expansion.pl
    ${projectDir}/bin/GO_make_isoform_hash.pl
    """
}
