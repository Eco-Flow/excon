process CONCAT_SINGLE_COPY {
    tag "$meta.id"
    label 'process_single'

    container "${ workflow.containerEngine == 'singularity' && !task.ext.singularity_pull_docker_container ?
        'https://depot.galaxyproject.org/singularity/eggnog-mapper:2.1.13--pyhdfd78af_2' :
        'biocontainers/eggnog-mapper:2.1.13--pyhdfd78af_2' }"

    input:
    tuple val(meta), path(orthofinder_dir)

    output:
    tuple val(meta), path('supermatrix.faa')   , emit: alignment
    tuple val(meta), path('partitions.txt')    , emit: partitions
    tuple val("${task.process}"), val('python'), eval("python3 --version | sed 's/Python //'"), emit: versions, topic: versions

    when:
    task.ext.when == null || task.ext.when

    script:
    def args = task.ext.args ?: ''
    """
    concat_single_copy.py \\
        -m ${orthofinder_dir}/MultipleSequenceAlignments \\
        -g ${orthofinder_dir}/Orthogroups/Orthogroups.tsv \\
        -o supermatrix.faa \\
        -p partitions.txt \\
        $args
    """

    stub:
    """
    touch supermatrix.faa
    touch partitions.txt
    """
}
