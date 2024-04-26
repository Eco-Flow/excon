process GFFREAD {
    label 'process_single'
    tag "$sample_id"
    container = 'ecoflowucl/gffread_python:python-3.11.9_Linux_x86_64_perl-5.36.0'

    input:
    tuple val(sample_id), path(fasta), path(gff)

    output:
    path( "${sample_id}.prot.fa" ), emit: proteins
    path( "${sample_id}.prot.fa.largestIsoform.fa" ), emit: longest
    path( "${sample_id}.splicedcds.fa" )
    path( "${sample_id}.splicedexons.fa" )
    path( "${sample_id}.gff_for_jvci.gff3" ), emit: gffs
    path( "${sample_id}_gene_alltran_list.txt" ), emit: gene_to_isoforms
    path( "${sample_id}.splicedcds.fa.nucl.longest.fa" )
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

    ${projectDir}/bin/gff_to_genetranshash.2.pl
    ${projectDir}/bin/prot_fasta_to_longest.pl ${sample_id}.prot.fa ${sample_id}_longestisoform.txt
    ${projectDir}/bin/fasta_topIsoform.pl ${sample_id}.splicedcds.fa ${sample_id}_longestisoform.txt


    #This part checks if longest isoform worked, if not we will continue with all proteins into Orthofinder. Warning sent to screen.
    #Largest isoforms has content if true
    #Largest isoforms does not have content if false. Just use full protein file (could be a genome without isoforms)

    if [[ -s ${sample_id}.prot.fa.largestIsoform.fa ]];then
      echo all_good
    else
      cp ${sample_id}.prot.fa ${sample_id}.prot.fa.largestIsoform.fa
    fi

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        gffread version: \$(gffread --version)
        Perl version: \$(perl --version | grep "version" | sed 's/.*(//g' | sed 's/[)].*//')
    END_VERSIONS
    """
}
