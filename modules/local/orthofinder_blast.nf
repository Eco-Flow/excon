process ORTHOFINDER_BLAST {
    tag "$meta.id"
    label 'process_high'
    label 'process_med_long'

    // The same OrthoFinder as the container. (modules/nf-core/orthofinder/environment.yml
    // pins 3.1.3, whose -op prints the search commands instead of writing them to a file.)
    conda "bioconda::orthofinder=3.2"
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
    # -op prepares the search databases and stops. OrthoFinder 3.2 writes the search
    # commands to WorkingDirectory/blast_commands.txt (3.1 printed them to stdout), one
    # complete command per line, for the program -S chose: DIAMOND, BLAST or MMseqs2.
    orthofinder \\
        -f input \\
        -n $prefix \\
        $args \\
        -op > orthofinder_prep.log

    # Run every search command, up to task.cpus at once. Each line goes into a script
    # of its own, so bash gets it unchanged (an MMseqs2 line is two commands joined
    # by ';').
    wd=input/OrthoFinder/Results_${prefix}/WorkingDirectory
    mkdir search_commands
    awk 'NF { f = sprintf("search_commands/%06d.sh", ++n); print > f; close(f) }' \$wd/blast_commands.txt
    n_cmds=\$(ls search_commands | wc -l)
    if [ "\$n_cmds" -eq 0 ]; then
        echo "ERROR: OrthoFinder -op wrote no search commands to \$wd/blast_commands.txt" >&2
        exit 1
    fi
    find search_commands -name '*.sh' | sort | xargs -P $task.cpus -n 1 bash

    # MMseqs2 leaves its databases (an index of about 0.9 GB per species, whatever
    # the proteome's size), a result database per pair and temporary directories.
    # The -b restart only needs the Blast*.txt results, so these are removed.
    rm -rf \$wd/mmseqsDB* \$wd/Blast*.txt.db* \$wd/tmpBlast*

    # The WorkingDirectory, now holding the Blast*.txt(.gz) results, is the output.
    # ORTHOFINDER_PHYLO restarts from it with -b; OrthoFinder 3.2 checks the searches
    # are all present and complete, and goes straight on to orthogroup inference.
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
