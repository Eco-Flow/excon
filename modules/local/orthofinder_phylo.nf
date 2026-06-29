process ORTHOFINDER_PHYLO {
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
    tuple val(meta), path(blast_wd)

    output:
    tuple val(meta), path("$prefix")                                            , emit: orthofinder
    tuple val(meta), path("$prefix/WorkingDirectory")                           , emit: working
    tuple val("${task.process}"), val('orthofinder'), eval("NO_COLOR=1 orthofinder --version | cut -d 'v' -f2 | perl -pe 's/\\e\\[[0-9;]*m//g'"), emit: versions_orthofinder, topic: versions
    path("$prefix/Orthogroups/Orthogroups.tsv")                                 , emit: orthologues
    path("$prefix/Species_Tree/SpeciesTree_rooted_node_labels.txt")             , emit: speciestree

    when:
    task.ext.when == null || task.ext.when

    script:
    def args = task.ext.args ?: ''
    prefix = task.ext.prefix ?: "${meta.id}"
    """
    # Build a local working directory from symlinks to the blast results.
    # This avoids any data copy (critical on Lustre/NFS where cp of many files is slow).
    # OrthoFinder writes new results (OrthoFinder/ subdir) into this local dir, so it
    # never touches the original blast work directory.
    # Any OrthoFinder/ dir from a previous run is explicitly excluded.
    mkdir local_blast_wd
    blast_abs=\$(readlink -f $blast_wd)
    for item in "\$blast_abs"/*; do
        name=\$(basename "\$item")
        [ "\$name" = "OrthoFinder" ] && continue
        ln -s "\$item" "local_blast_wd/\$name"
    done

    orthofinder \\
        -t $task.cpus \\
        -a ${[task.cpus, 4].min()} \\
        -b local_blast_wd \\
        -n $prefix \\
        $args

    # OrthoFinder writes results to local_blast_wd/OrthoFinder/Results_$prefix
    mv local_blast_wd/OrthoFinder/Results_${prefix} $prefix
    """

    stub:
    prefix = task.ext.prefix ?: "${meta.id}"
    """
    mkdir -p    $prefix/Comparative_Genomics_Statistics
    mkdir       $prefix/Gene_Duplication_Events
    mkdir       $prefix/Gene_Trees
    mkdir       $prefix/Orthogroup_Sequences
    mkdir       $prefix/Orthogroups
    mkdir       $prefix/Orthologues
    mkdir       $prefix/Phylogenetic_Hierarchical_Orthogroups
    mkdir       $prefix/Phylogenetically_Misplaced_Genes
    mkdir       $prefix/Putative_Xenologs
    mkdir       $prefix/Resolved_Gene_Trees
    mkdir       $prefix/Single_Copy_Orthologue_Sequences
    mkdir       $prefix/Species_Tree
    mkdir       $prefix/WorkingDirectory
    touch       $prefix/Log.txt
    touch       $prefix/Orthogroups/Orthogroups.tsv
    touch       $prefix/Species_Tree/SpeciesTree_rooted_node_labels.txt
    """
}
