process EXTRACT_SINGLE_COPY {
    tag "$meta.id"
    label 'process_single'

    container "${ workflow.containerEngine == 'singularity' && !task.ext.singularity_pull_docker_container ?
        'https://depot.galaxyproject.org/singularity/eggnog-mapper:2.1.13--pyhdfd78af_2' :
        'biocontainers/eggnog-mapper:2.1.13--pyhdfd78af_2' }"

    input:
    tuple val(meta), path(orthofinder_dir), path(proteomes, stageAs: 'proteomes/*')

    output:
    tuple val(meta), path('single_copy_orthogroups'), emit: orthogroups
    tuple val("${task.process}"), val('python'), eval("python3 --version | sed 's/Python //'"), emit: versions, topic: versions

    when:
    task.ext.when == null || task.ext.when

    script:
    def args = task.ext.args ?: ''
    """
    extract_single_copy.py \\
        -g ${orthofinder_dir}/Orthogroups/Orthogroups.tsv \\
        -p proteomes \\
        -o single_copy_orthogroups \\
        $args
    """

    stub:
    """
    mkdir single_copy_orthogroups
    """
}
