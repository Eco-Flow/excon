process GFFREAD {
    label 'process_single'
    tag "$sample_id"
    container = 'ecoflowucl/gffread_python:python-3.11.9_Linux_x86_64_perl-5.36.0'

    input:
    tuple val(sample_id), path(fasta), path(gff)

    output:
    path( "${sample_id}.prot.fa" ), emit: proteins
    tuple val(sample_id), path("${sample_id}.prot.fa" ), emit: proteins_busco
    path( "${sample_id}.prot.fa" ), emit: longest
    path( "${sample_id}.splicedcds.fa" )
    path( "${sample_id}.splicedexons.fa" )
    path( "${sample_id}.gff_for_jvci.gff3" ), emit: gffs
    tuple val(sample_id), path("${sample_id}.gff_for_jvci.gff3"), emit: gffs_agat
    path( "${sample_id}_gene_alltran_list.txt" ), emit: gene_to_isoforms
    path( "${sample_id}.splicedcds.fa" )
    tuple val( "${sample_id}" ), path( "${fasta}" ), emit: fasta_quast
    path "versions.yml", emit: versions

    script:
    """
    #Check if gff3 or genome file is gzipped:
    if [[ $gff == *.gz ]]; then
      zcat $gff > gff_temp
    else
      cp $gff gff_temp
    fi

    if [[ $fasta == *.gz ]]; then
      zcat $fasta > genome_temp
    else
      cp $fasta genome_temp
    fi

    #Convert Augustus gff files if found, then do gffread to print out the nucleotide files for each gene.
    # Note: The GFF input now comes from LONGEST process (agat longest isoform), so it should already be processed

    head -n 1 gff_temp > tbd

    if grep -q AUGUSTUS tbd; then 
        gtf2gff.pl <gff_temp --out=${sample_id}.gff_for_jvci.gff3
        #Remove lines of the GFF that have ? in the strand section, as this cannot be parsed by gffread
        awk '\$7 != "?" { print \$0 }' ${sample_id}.gff_for_jvci.gff3  > ${sample_id}.gff_for_jvci.noquest.gff3
        
        gffread -w ${sample_id}.splicedexons.fa -g genome_temp ${sample_id}.gff_for_jvci.noquest.gff3
        gffread -x ${sample_id}.splicedcds.fa -g genome_temp ${sample_id}.gff_for_jvci.noquest.gff3
        gffread -y ${sample_id}.prot.fa -g genome_temp ${sample_id}.gff_for_jvci.noquest.gff3 -F -S

    else
        mv gff_temp ${sample_id}.gff_for_jvci.gff3
        #Remove lines of the GFF that have ? in the strand section, as this cannot be parsed by gffread
        awk '\$7 != "?" { print \$0 }' ${sample_id}.gff_for_jvci.gff3  > ${sample_id}.gff_for_jvci.noquest.gff3
        
        gffread -w ${sample_id}.splicedexons.fa -g genome_temp ${sample_id}.gff_for_jvci.noquest.gff3
        gffread -x ${sample_id}.splicedcds.fa -g genome_temp ${sample_id}.gff_for_jvci.noquest.gff3
        gffread -y ${sample_id}.prot.fa -g genome_temp ${sample_id}.gff_for_jvci.noquest.gff3 -F -S

    fi

    # Create gene to isoform mapping (still needed for downstream processes)
    ${projectDir}/bin/gff_to_genetranshash.2.pl

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        gffread version: \$(gffread --version)
        Perl version: \$(perl --version | grep "version" | sed 's/.*(//g' | sed 's/[)].*//')
    END_VERSIONS
    """
}
