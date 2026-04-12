process NCBIGENOMEDOWNLOAD {
    tag "$meta.id"
    label 'process_single'

    conda "${moduleDir}/environment.yml"
    container "${ workflow.containerEngine == 'singularity' && !task.ext.singularity_pull_docker_container ?
        'https://depot.galaxyproject.org/singularity/ncbi-genome-download:0.3.3--pyh7cba7a3_0' :
        'biocontainers/ncbi-genome-download:0.3.3--pyh7cba7a3_0' }"

    input:
    val meta
    path accessions
    path taxids
    val groups

    output:
    tuple val(meta), path("*_genomic.gbff.gz"),         emit: gbk,      optional: true
    tuple val(meta), path("*_genomic.fna.gz"),          emit: fna,      optional: true
    tuple val(meta), path("*_rm.out.gz"),               emit: rm,       optional: true
    tuple val(meta), path("*_feature_table.txt.gz"),    emit: features, optional: true
    tuple val(meta), path("*_genomic.gff.gz"),          emit: gff,      optional: true
    tuple val(meta), path("*_protein.faa.gz"),          emit: faa,      optional: true
    tuple val(meta), path("*_protein.gpff.gz"),         emit: gpff,     optional: true
    tuple val(meta), path("*_wgsmaster.gbff.gz"),       emit: wgs_gbk,  optional: true
    tuple val(meta), path("*_cds_from_genomic.fna.gz"), emit: cds,      optional: true
    tuple val(meta), path("*_rna.fna.gz"),              emit: rna,      optional: true
    tuple val(meta), path("*_rna_from_genomic.fna.gz"), emit: rna_fna,  optional: true
    tuple val(meta), path("*_assembly_report.txt"),     emit: report,   optional: true
    tuple val(meta), path("*_assembly_stats.txt"),      emit: stats,    optional: true
    tuple val("${task.process}"), val('ncbigenomedownload'), eval('ncbi-genome-download --version'), topic: versions, emit: versions_ncbigenomedownload

    when:
    task.ext.when == null || task.ext.when

    script:
def args           = task.ext.args ?: ''
def accessions_opt = accessions ? "-A accessions.clean.txt" : ""
def taxids_opt     = taxids ? "-t ${taxids}" : ""

"""
set -euo pipefail
shopt -s nullglob

if [ -f "${accessions}" ]; then
    awk '{\$1=\$1; print}' "${accessions}" > accessions.clean.txt
    echo "Cleaned accessions:"
    cat accessions.clean.txt
fi

echo "NGD version:"
ncbi-genome-download --version

ncbi-genome-download \\
    $args \\
    $accessions_opt \\
    $taxids_opt \\
    --output-folder ./ \\
    --flat-output \\
    --parallel 1 \\
    $groups

gz_files=( *.gz )

if [ \${#gz_files[@]} -eq 0 ]; then
    echo "ERROR: ncbi-genome-download produced no .gz files" >&2
    ls -lah >&2 || true
    exit 1
fi

for f in "\${gz_files[@]}"; do
    if ! gzip -t "\$f" 2>/dev/null; then
        echo "ERROR: '\$f' is not valid gzip" >&2
        head -c 500 "\$f" >&2 || true
        exit 1
    fi
done
"""



}

