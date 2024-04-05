process GFFREAD {
    label 'gffread'
    tag "$sample_id"
    container = 'chriswyatt/bioseqio_gffread'
    publishDir "$params.outdir/Gffread" , mode: "copy"
             
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

    script:
    """
    #Check if gff3 or genome file is gzipped:
if [[ $gff == *.gz ]]
then
    zcat $gff > gff_temp
else
    cp $gff gff_temp
fi


if [[ $fasta == *.gz ]]
then
    zcat $fasta > genome_temp
else
    cp $fasta genome_temp
fi

    #Convert Augustus gff files if found, then do gffread to print out the nucleotide files for each gene.

    head -n 1 gff_temp > tbd

    if grep -q AUGUSTUS tbd; then 
        #python3 $projectDir/bin/convert_augustus_to_gffs.py -i gff_temp -o ${sample_id}.gff_for_jvci.gff3
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

	gff_to_genetranshash.2.pl
	prot_fasta_to_longest.pl ${sample_id}.prot.fa ${sample_id}_longestisoform.txt
	fasta_topIsoform.pl ${sample_id}.splicedcds.fa ${sample_id}_longestisoform.txt	     
    """
}


    
