process QUAST {
    //tag "$meta.id"
    label 'process_medium'

    conda "${moduleDir}/environment.yml"
    container "${ workflow.containerEngine == 'singularity' && !task.ext.singularity_pull_docker_container ?
        'https://depot.galaxyproject.org/singularity/quast:5.2.0--py39pl5321h2add14b_1' :
        'biocontainers/quast:5.2.0--py39pl5321h2add14b_1' }"

    input:
    tuple val(meta) , path(consensus)
    //tuple val(meta2), path(fasta)
    tuple val(meta3), path(gff)

    output:
    tuple val(meta), path("${meta}/*")                 , emit: results
    tuple val(meta), path("${meta}.tsv")               , emit: tsv
    tuple val(meta), path("${meta}_transcriptome.tsv") , optional: true , emit: transcriptome
    tuple val(meta), path("${meta}_misassemblies.tsv") , optional: true , emit: misassemblies
    tuple val(meta), path("${meta}_unaligned.tsv")     , optional: true , emit: unaligned
    path "versions.yml"                                , emit: versions

    when:
    task.ext.when == null || task.ext.when

    script:
    def args      = task.ext.args   ?: ''
    //prefix        = task.ext.prefix ?: "${meta.id}"
    def features  = gff             ?  "--features $gff" : ''
    //def reference = fasta           ?  "-r $fasta"       : ''
    """
    quast.py \\
        --output-dir $meta \\
        $features \\
        --threads $task.cpus \\
        $args \\
        ${consensus.join(' ')}

    ln -s ${meta}/report.tsv ${meta}.tsv
    [ -f  ${meta}/contigs_reports/all_alignments_transcriptome.tsv ] && ln -s ${meta}/contigs_reports/all_alignments_transcriptome.tsv ${meta}_transcriptome.tsv
    [ -f  ${meta}/contigs_reports/misassemblies_report.tsv         ] && ln -s ${meta}/contigs_reports/misassemblies_report.tsv ${meta}_misassemblies.tsv
    [ -f  ${meta}/contigs_reports/unaligned_report.tsv             ] && ln -s ${meta}/contigs_reports/unaligned_report.tsv ${meta}_unaligned.tsv

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        quast: \$(quast.py --version 2>&1 | sed 's/^.*QUAST v//; s/ .*\$//')
    END_VERSIONS
    """

    stub:
    def args      = task.ext.args   ?: ''
    //prefix        = task.ext.prefix ?: "${meta.id}"
    def features  = gff             ? "--features $gff" : ''
    //def reference = fasta           ? "-r $fasta" : ''

    """
    mkdir -p $meta
    touch $meta/report.tsv
    touch $meta/report.html
    touch $meta/report.pdf
    touch $meta/quast.log
    touch $meta/transposed_report.txt
    touch $meta/transposed_report.tex
    touch $meta/icarus.html
    touch $meta/report.tex
    touch $meta/report.txt

    mkdir -p $meta/basic_stats
    touch $meta/basic_stats/cumulative_plot.pdf
    touch $meta/basic_stats/Nx_plot.pdf
    touch $meta/basic_stats/genome_GC_content_plot.pdf
    touch $meta/basic_stats/GC_content_plot.pdf

    mkdir -p $meta/icarus_viewers
    touch $meta/icarus_viewers/contig_size_viewer.html

    ln -s $meta/report.tsv ${meta}.tsv

    if [ $fasta ]; then
        touch $meta/basic_stats/NGx_plot.pdf
        touch $meta/basic_stats/gc.icarus.txt

        mkdir -p $meta/aligned_stats
        touch $meta/aligned_stats/NAx_plot.pdf
        touch $meta/aligned_stats/NGAx_plot.pdf
        touch $meta/aligned_stats/cumulative_plot.pdf

        mkdir -p $meta/contigs_reports
        touch $meta/contigs_reports/all_alignments_transcriptome.tsv
        touch $meta/contigs_reports/contigs_report_transcriptome.mis_contigs.info
        touch $meta/contigs_reports/contigs_report_transcriptome.stderr
        touch $meta/contigs_reports/contigs_report_transcriptome.stdout
        touch $meta/contigs_reports/contigs_report_transcriptome.unaligned.info
        mkdir -p $meta/contigs_reports/minimap_output
        touch $meta/contigs_reports/minimap_output/transcriptome.coords
        touch $meta/contigs_reports/minimap_output/transcriptome.coords.filtered
        touch $meta/contigs_reports/minimap_output/transcriptome.coords_tmp
        touch $meta/contigs_reports/minimap_output/transcriptome.sf
        touch $meta/contigs_reports/minimap_output/transcriptome.unaligned
        touch $meta/contigs_reports/minimap_output/transcriptome.used_snps
        touch $meta/contigs_reports/misassemblies_frcurve_plot.pdf
        touch $meta/contigs_reports/misassemblies_plot.pdf
        touch $meta/contigs_reports/misassemblies_report.tex
        touch $meta/contigs_reports/misassemblies_report.tsv
        touch $meta/contigs_reports/misassemblies_report.txt
        touch $meta/contigs_reports/transcriptome.mis_contigs.fa
        touch $meta/contigs_reports/transposed_report_misassemblies.tex
        touch $meta/contigs_reports/transposed_report_misassemblies.tsv
        touch $meta/contigs_reports/transposed_report_misassemblies.txt
        touch $meta/contigs_reports/unaligned_report.tex
        touch $meta/contigs_reports/unaligned_report.tsv
        touch $meta/contigs_reports/unaligned_report.txt

        mkdir -p $meta/genome_stats
        touch $meta/genome_stats/genome_info.txt
        touch $meta/genome_stats/transcriptome_gaps.txt
        touch $meta/icarus_viewers/alignment_viewer.html

        ln -sf ${meta}/contigs_reports/misassemblies_report.tsv ${meta}_misassemblies.tsv
        ln -sf ${meta}/contigs_reports/unaligned_report.tsv ${meta}_unaligned.tsv
        ln -sf ${meta}/contigs_reports/all_alignments_transcriptome.tsv ${meta}_transcriptome.tsv

    fi

    if ([ $fasta ] && [ $gff ]); then
        touch $meta/genome_stats/features_cumulative_plot.pdf
        touch $meta/genome_stats/features_frcurve_plot.pdf
    fi

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        quast: \$(quast.py --version 2>&1 | sed 's/^.*QUAST v//; s/ .*\$//')
    END_VERSIONS
    """
}
