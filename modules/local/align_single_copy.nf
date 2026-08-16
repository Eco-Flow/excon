process ALIGN_SINGLE_COPY {
    tag "$meta.id"
    label 'process_medium'

    conda "${moduleDir}/../nf-core/mafft/align/environment.yml"
    container "${ workflow.containerEngine == 'singularity' && !task.ext.singularity_pull_docker_container ?
        'https://depot.galaxyproject.org/singularity/mulled-v2-12eba4a074f913c639117640936668f5a6a01da6:425707898cf4f85051b77848be253b88f1d2298a-0' :
        'quay.io/biocontainers/mulled-v2-12eba4a074f913c639117640936668f5a6a01da6:425707898cf4f85051b77848be253b88f1d2298a-0' }"

    input:
    tuple val(meta), path(orthogroup_dir)

    output:
    tuple val(meta), path('single_copy_alignments'), emit: alignments
    tuple val("${task.process}"), val('mafft'), eval("mafft --version 2>&1 | sed 's/^v//;s/ .*//'"), emit: versions, topic: versions

    when:
    task.ext.when == null || task.ext.when

    script:
    def args = task.ext.args ?: '--auto'
    // One task aligning every orthogroup internally, rather than one Nextflow task
    // per orthogroup: the orthogroup count scales with the number of species and
    // would otherwise swamp the scheduler.
    """
    mkdir single_copy_alignments

    # -L because Nextflow stages the input directory as a symlink, which
    # plain find will not descend into.
    find -L ${orthogroup_dir} -name '*.fa' | sort > single_copy_files.txt

    if [ ! -s single_copy_files.txt ]; then
        echo "ERROR: no single-copy orthogroup sequences found in ${orthogroup_dir}" >&2
        exit 1
    fi

    xargs -a single_copy_files.txt -P ${task.cpus} -I {} \\
        sh -c 'mafft $args --quiet --anysymbol "\$1" > "single_copy_alignments/\$(basename "\$1")"' _ {}

    echo "Aligned \$(ls single_copy_alignments | wc -l) of \$(wc -l < single_copy_files.txt) single-copy orthogroups"
    """

    stub:
    """
    mkdir single_copy_alignments
    """
}
