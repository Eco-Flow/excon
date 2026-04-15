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
    # Prepare DIAMOND databases and blast_commands.txt (-op stops before running searches)
    orthofinder \\
        -f input \\
        -n $prefix \\
        $args \\
        -op

    mv input/OrthoFinder/Results_${prefix}/WorkingDirectory ${meta.id}_blast_wd

    # Run the DIAMOND searches prepared by OrthoFinder
    # Each command writes its Blast*.txt.gz output into the WorkingDirectory
    ( cd ${meta.id}_blast_wd && bash blast_commands.txt )

    # Remove blast_commands.txt so orthofinder -b knows searches are complete
    # (orthofinder -b only re-runs searches if this file is present)
    rm ${meta.id}_blast_wd/blast_commands.txt
    """

    stub:
    prefix = task.ext.prefix ?: "${meta.id}"
    """
    mkdir -p ${meta.id}_blast_wd
    touch ${meta.id}_blast_wd/SpeciesIDs.txt
    touch ${meta.id}_blast_wd/SequenceIDs.txt
    """
}
