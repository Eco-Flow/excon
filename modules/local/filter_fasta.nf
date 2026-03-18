process FILTER_FASTA {
    tag "$meta.id"
    label 'process_single'

    input:
    tuple val(meta), path(fasta)

    output:
    tuple val(meta), path("${meta.id}.clean.fasta")

    script:
    """
    sed '/^[^>]/s/[.*]//g' ${fasta} > ${meta.id}.clean.fasta
    """
}