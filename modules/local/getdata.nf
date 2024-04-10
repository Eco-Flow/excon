process GET_DATA {
    label 'process_low'

    container = 'ecoflowucl/biomart_perl:r-4.3.3_perl-5.38.2'
    
    input:
    val(ensembl_repo)
    val(ensembl_host)
    val(public_species)
        
    output:
    path("${public_species}.fasta") , emit: fasta_files
    path("${public_species}.go.txt") , emit: gene_ontology_files
    path "versions.yml", emit: versions

    script:
    """
    #Pull all Biomart records for species.
    ${projectDir}/bin/R_biomart.R '$ensembl_repo' '$ensembl_host' '$public_species'
    #Tidy up records
    #Insert custom perl/unix script.
    cat Myoutput* > All_fasta
    sed '/ensembl_gene_id/d' All_fasta > All_fasta2
    sed '/peptide/d' All_fasta2 > All_fasta3
    ${projectDir}/bin/Fix_fasta.pl All_fasta3 > ${public_species}.fasta
    mv go_hash.txt ${public_species}.go.txt

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        R version: \$(R --version | grep "R version" | sed 's/[(].*//' | sed 's/ //g' | sed 's/[^0-9]*//')
        Perl version: \$(perl --version | grep "version" | sed 's/.*(//g' | sed 's/[)].*//')
    END_VERSIONS
    """
}
