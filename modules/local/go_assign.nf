process GO_ASSIGN {
    label 'process_single'
    label 'error_ignore'
    tag "GO_assign_${Focal}"
    container = 'ecoflowucl/chopgo:r-4.3.2_python-3.10_perl-5.38'

    input:
    path Go_files
    path Orthogroups
    path Focal
    path Gene_to_trans

    output:
    path("${Focal}_Sort_Result_All_Combine_GO_format") , emit: go_hash
    path("OG_GO_format.tsv") , emit: go_og
    path("*results_ALL.tab.pdf") , emit: duplicate_go, optional: true
    path("*family_expansions.txt") , emit: go_counts
    path("*transcripts_Combine_GO_format.txt"), emit: trans_go
    path "versions.yml", emit: versions

    script:
    """
    ${projectDir}/bin/Goatee_ortho_go_match.pl ${Orthogroups} ${Focal}
    ${projectDir}/bin/Orthofinder_duplicate_go.pl
    ${projectDir}/bin/Orthofinder_gene_expansion.pl
    ${projectDir}/bin/GO_make_isoform_hash.pl

    #Sort output of the GO assignment, so we have consistent order.
    sort ${Focal}_Result_All_Combine_GO_format > ${Focal}_Sort_Result_All_Combine_GO_format

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        R version: \$(R --version | grep "R version" | sed 's/[(].*//' | sed 's/ //g' | sed 's/[^0-9]*//')
        GO version: \$(Rscript -e "as.data.frame(installed.packages())[ ,c(1,3)]" | grep topGO | sed 's/[^0-9]*//')
        Perl version: \$(perl --version | grep "version" | sed 's/.*(//g' | sed 's/[)].*//')
    END_VERSIONS
    """
}
