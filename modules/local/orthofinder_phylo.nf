process ORTHOFINDER_PHYLO {
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
    # Build the working directory from symlinks to the search results, at
    # $prefix/WorkingDirectory: no data is copied (critical on Lustre/NFS), and the
    # BLAST task's own directory is never written to. On a -b restart OrthoFinder 3.2
    # writes its results into the parent of the directory it is given, so they land
    # in $prefix/ (3.1 created WorkingDirectory/OrthoFinder/Results_<name>/ instead).
    # Any OrthoFinder/ dir from a previous run is left out. Subdirectories (3.2's
    # dependencies/, a few KB) are copied: OrthoFinder rewrites files in them, which
    # through a symlink would change the BLAST task's own output.
    mkdir -p $prefix/WorkingDirectory
    blast_abs=\$(readlink -f $blast_wd)
    for item in "\$blast_abs"/*; do
        name=\$(basename "\$item")
        [ "\$name" = "OrthoFinder" ] && continue
        if [ -d "\$item" ]; then
            cp -rL "\$item" "$prefix/WorkingDirectory/\$name"
        else
            ln -s "\$item" "$prefix/WorkingDirectory/\$name"
        fi
    done

    # OrthoFinder 3's "Initial processing of each species" step has a 200 s
    # STALL_TIMEOUT that misfires as a failure on any species whose own processing
    # genuinely takes longer than that -- not an actual hang. No CLI flag exists to
    # raise it (davidemms/OrthoFinder#1024), and the container's site-packages is
    # read-only, so patch the function's default via a sitecustomize.py on PYTHONPATH.
    # It is the last default argument of
    #   3.2:   orthofinder.mcl_orthogroups.infer_orthogroups.do_orthogroups
    #   3.1.x: orthofinder.orthogroups.gathering.DoOrthogroups
    # The `orthofinder` package isn't on the normal Python path -- the launcher adds
    # its own src/ dir at run time, after Python's startup already tried to import
    # sitecustomize.py -- so add that src/ dir to PYTHONPATH ourselves.
    orthofinder_src="\$(dirname "\$(command -v orthofinder)")/src"
    mkdir pypatch
    cat > pypatch/sitecustomize.py <<'PYEOF'
try:
    from orthofinder.mcl_orthogroups import infer_orthogroups as _mod
    _fn = _mod.do_orthogroups
except ImportError:
    from orthofinder.orthogroups import gathering as _mod
    _fn = _mod.DoOrthogroups
_fn.__defaults__ = _fn.__defaults__[:-1] + (3600.,)
PYEOF
    export PYTHONPATH="\${orthofinder_src}:\$PWD/pypatch\${PYTHONPATH:+:\$PYTHONPATH}"

    # Fail fast if the patch didn't take: sitecustomize.py import errors are swallowed
    # by Python's startup (a warning, nothing more), so check the default really is 3600 s
    python3 -c "
import sitecustomize
try:
    from orthofinder.mcl_orthogroups.infer_orthogroups import do_orthogroups as f
except ImportError:
    from orthofinder.orthogroups.gathering import DoOrthogroups as f
assert f.__defaults__[-1] == 3600., f.__defaults__
" || {
        echo "ERROR: the OrthoFinder STALL_TIMEOUT patch did not apply ('\${orthofinder_src}'). The container's package layout may have changed." >&2
        exit 1
    }

    orthofinder \\
        -t $task.cpus \\
        -a ${[task.cpus, 4].min()} \\
        -b $prefix/WorkingDirectory \\
        -n $prefix \\
        $args

    if [ ! -s $prefix/Orthogroups/Orthogroups.tsv ]; then
        echo "ERROR: OrthoFinder wrote no $prefix/Orthogroups/Orthogroups.tsv; see its output above." >&2
        exit 1
    fi
    """

    stub:
    prefix = task.ext.prefix ?: "${meta.id}"
    """
    mkdir -p    $prefix/Comparative_Genomics_Statistics
    mkdir       $prefix/Gene_Duplication_Events
    mkdir       $prefix/Gene_Trees
    mkdir       $prefix/Orthogroup_Sequences
    mkdir       $prefix/Orthogroups
    mkdir       $prefix/Orthologs
    mkdir       $prefix/Phylogenetic_Hierarchical_Orthogroups
    mkdir       $prefix/Phylogenetically_Misplaced_Genes
    mkdir       $prefix/Putative_Xenologs
    mkdir       $prefix/Resolved_Gene_Trees
    mkdir       $prefix/Single_Copy_Ortholog_Sequences
    mkdir       $prefix/MultipleSequenceAlignments
    mkdir       $prefix/Species_Tree
    mkdir       $prefix/WorkingDirectory
    touch       $prefix/Log.txt
    touch       $prefix/Orthogroups/Orthogroups.tsv
    touch       $prefix/Species_Tree/SpeciesTree_rooted_node_labels.txt
    """
}
