process ORTHOFINDER_ASSIGN {
    tag "$meta.id"
    label 'process_high'
    label 'process_med_long'

    conda "bioconda::orthofinder=3.2"
    container {
        workflow.containerEngine == 'singularity' && !task.ext?.singularity_pull_docker_container ?
        'https://depot.galaxyproject.org/singularity/orthofinder:3.2--hdfd78af_0' :
        'biocontainers/orthofinder:3.2--hdfd78af_0'
    }

    input:
    // core: the results of the OrthoFinder run on the core species (ORTHOFINDER_PHYLO)
    // fastas: the proteomes of every other species, to assign to its orthogroups
    tuple val(meta), path(core)
    path(fastas, stageAs: 'assign/*')

    output:
    tuple val(meta), path("$prefix")                                            , emit: orthofinder
    tuple val("${task.process}"), val('orthofinder'), eval("NO_COLOR=1 orthofinder --version | cut -d 'v' -f2 | perl -pe 's/\\e\\[[0-9;]*m//g'"), emit: versions_orthofinder, topic: versions
    path("$prefix/Orthogroups/Orthogroups.tsv")                                 , emit: orthologues
    path("$prefix/Species_Tree/SpeciesTree_rooted_node_labels.txt")             , emit: speciestree

    when:
    task.ext.when == null || task.ext.when

    script:
    def args = task.ext.args ?: ''
    prefix = task.ext.prefix ?: "${meta.id}_assign"
    def n_assign = fastas instanceof List ? fastas.size() : 1
    """
    # The same STALL_TIMEOUT patch as ORTHOFINDER_PHYLO: --assign clusters the genes it
    # can't place in a core orthogroup with the standard workflow, whose "Initial
    # processing" step has OrthoFinder 3's 200 s timer (davidemms/OrthoFinder#1024).
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

    # --assign adds each species in assign/ to the core run's orthogroups (SHOOT
    # profiles, no all-versus-all search) and infers the species tree with ASTRAL-Pro.
    # It reads the core results without changing them, and writes Results_<-n> into
    # the parent of --core: here the task directory, since $core is a link to them.
    orthofinder \\
        --assign assign \\
        --core $core \\
        -t $task.cpus \\
        -a ${[task.cpus, 4].min()} \\
        -n $prefix \\
        $args

    mv Results_$prefix $prefix
    if [ ! -s $prefix/Orthogroups/Orthogroups.tsv ]; then
        echo "ERROR: OrthoFinder wrote no $prefix/Orthogroups/Orthogroups.tsv; see its output above." >&2
        exit 1
    fi

    # CAFE needs every species: check the core's and the assigned ones are all there
    n_core=\$(head -1 $core/Orthogroups/Orthogroups.tsv | awk -F'\\t' '{ print NF - 1 }')
    n_all=\$(head -1 $prefix/Orthogroups/Orthogroups.tsv | awk -F'\\t' '{ print NF - 1 }')
    if [ "\$n_all" -ne \$(( n_core + $n_assign )) ]; then
        echo "ERROR: $prefix/Orthogroups/Orthogroups.tsv has \$n_all species, expected \$n_core core + $n_assign assigned." >&2
        exit 1
    fi
    """

    stub:
    prefix = task.ext.prefix ?: "${meta.id}_assign"
    """
    mkdir -p $prefix/Orthogroups $prefix/Species_Tree $prefix/WorkingDirectory
    touch    $prefix/Orthogroups/Orthogroups.tsv
    touch    $prefix/Species_Tree/SpeciesTree_rooted_node_labels.txt
    """
}
