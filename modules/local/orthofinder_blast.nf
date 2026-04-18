process ORTHOFINDER_BLAST {
    tag "$meta.id"
    label 'process_high'
    label 'process_med_long'

    conda "${moduleDir}/../nf-core/orthofinder/environment.yml"
    container {
        workflow.containerEngine == 'singularity' && !task.ext?.singularity_pull_docker_container ?
        'https://depot.galaxyproject.org/singularity/orthofinder:3.1.4--hdfd78af_0' :
        'biocontainers/orthofinder:3.1.4--hdfd78af_0'
    }

    input:
    tuple val(meta), path(fastas, stageAs: 'input/')

    output:
    tuple val(meta), path("${meta.id}_blast_wd"), emit: working_dir
    tuple val("${task.process}"), val('orthofinder'), eval("NO_COLOR=1 orthofinder --version | cut -d 'v' -f2 | perl -pe 's/\\e\\[[0-9;]*m//g'"), emit: versions_orthofinder, topic: versions

    when:
    task.ext.when == null || task.ext.when

    script:
    def args = task.ext.args ?: ''
    prefix = task.ext.prefix ?: "${meta.id}"
    """
    # -op prepares DIAMOND databases and prints commands to stdout (v3 does not write
    # blast_commands.txt — commands appear on stdout after "diamond commands that must be run")
    orthofinder \\
        -f input \\
        -n $prefix \\
        $args \\
        -op > orthofinder_prep.log

    # Run diamond commands in parallel (up to task.cpus at once).
    # Paths are relative to the current work directory.
    grep "^diamond blastp" orthofinder_prep.log | xargs -P $task.cpus -I {} bash -c '{}'

    # Move the WorkingDirectory (now containing Blast*.txt.gz) to the output location.
    # Because blast_commands.txt was never written, orthofinder -b will skip re-running
    # searches and proceed directly to orthogroup inference.
    mv input/OrthoFinder/Results_${prefix}/WorkingDirectory ${meta.id}_blast_wd
    """

    stub:
    prefix = task.ext.prefix ?: "${meta.id}"
    """
    mkdir -p ${meta.id}_blast_wd
    touch ${meta.id}_blast_wd/SpeciesIDs.txt
    touch ${meta.id}_blast_wd/SequenceIDs.txt
    """
}
