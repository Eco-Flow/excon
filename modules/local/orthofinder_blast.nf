process ORTHOFINDER_BLAST {
    tag "$meta.id"
    label 'process_high'
    label 'process_med_long'

    conda "${moduleDir}/../nf-core/orthofinder/environment.yml"
    container {
        workflow.containerEngine == 'singularity' && !task.ext?.singularity_pull_docker_container ?
        'https://depot.galaxyproject.org/singularity/orthofinder:3.2--hdfd78af_0' :
        'biocontainers/orthofinder:3.2--hdfd78af_0'
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
    # -op prepares the DIAMOND databases and stops. OrthoFinder 3.2 writes the search
    # commands to WorkingDirectory/blast_commands.txt (3.1 printed them to stdout).
    orthofinder \\
        -f input \\
        -n $prefix \\
        $args \\
        -op > orthofinder_prep.log

    wd=input/OrthoFinder/Results_${prefix}/WorkingDirectory
    n_cmds=\$(grep -c "^diamond blastp" \$wd/blast_commands.txt || true)
    if [ "\$n_cmds" -eq 0 ]; then
        echo "ERROR: OrthoFinder -op wrote no DIAMOND commands to \$wd/blast_commands.txt" >&2
        exit 1
    fi

    # Run the DIAMOND searches in parallel, up to task.cpus at once
    grep "^diamond blastp" \$wd/blast_commands.txt | xargs -P $task.cpus -I {} bash -c '{}'

    # The WorkingDirectory, now holding Blast*.txt.gz, is the output. ORTHOFINDER_PHYLO
    # restarts from it with -b; OrthoFinder 3.2 checks the searches are all present
    # and complete, and goes straight on to orthogroup inference.
    mv \$wd ${meta.id}_blast_wd
    """

    stub:
    prefix = task.ext.prefix ?: "${meta.id}"
    """
    mkdir -p ${meta.id}_blast_wd
    touch ${meta.id}_blast_wd/SpeciesIDs.txt
    touch ${meta.id}_blast_wd/SequenceIDs.txt
    """
}
