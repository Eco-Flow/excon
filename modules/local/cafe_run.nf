process CAFE_RUN {
    tag "${meta.id}"
    label 'process_medium'
    label 'process_long'
    container 'ecoflowucl/cafe:r-4.3.1'

    input:
    tuple val(meta), path(counts), path(tree)

    output:
    tuple val(meta), path("Out_${meta.id}"), emit: results
    tuple val("${task.process}"), val('cafe'), val('4.2.1'), emit: versions_cafe, topic: versions

    script:
    def args = task.ext.args ?: ''
    """
    cafe5 -i ${counts} -t ${tree} --cores ${task.cpus} -o Out_${meta.id} ${args}
    """
}
