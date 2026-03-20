process FILTER_FASTA {
    tag "$meta.id"
    label 'process_single'

    conda "conda-forge::sed=4.8"
    container "${ workflow.containerEngine == 'singularity' && !task.ext.singularity_pull_docker_container ?
        'https://community-cr-prod.seqera.io/docker/registry/v2/blobs/sha256/52/52ccce28d2ab928ab862e25aae26314d69c8e38bd41ca9431c67ef05221348aa/data' :
        'community.wave.seqera.io/library/coreutils_grep_gzip_lbzip2_pruned:838ba80435a629f8' }"


    input:
    tuple val(meta), path(fasta)

    output:
    tuple val(meta), path("${meta.id}.clean.fasta")

    script:
    """
    sed '/^[^>]/s/[.*]//g' ${fasta} > ${meta.id}.clean.fasta
    """
}
