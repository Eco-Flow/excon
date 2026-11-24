process PREPARE_PROTEINS {
    tag "$meta.id"
    label 'process_single'

    container "${ workflow.containerEngine == 'singularity' && !task.ext.singularity_pull_docker_container ?
        'https://depot.galaxyproject.org/singularity/eggnog-mapper:2.1.13--pyhdfd78af_2' :
        'biocontainers/eggnog-mapper:2.1.13--pyhdfd78af_2' }"

    input:
    tuple val(meta), path(proteins)

    output:
    tuple val(meta), path("${meta.id}.clean.fasta"),              emit: fasta
    tuple val(meta), path("${meta.id}.protein_ids.gff"),          emit: gff
    tuple val(meta), path("${meta.id}.internal_stop_codons.tsv"), emit: internal_stop_report, optional: true
    tuple val("${task.process}"), val('python'), eval("python3 --version | sed 's/Python //'"), emit: versions, topic: versions

    script:
    """
    python3 ${projectDir}/bin/prepare_proteins.py ${meta.id} ${proteins} ${params.internal_stop_action}
    """
}
