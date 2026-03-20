process EGGNOG_DOWNLOAD {
    tag "eggnog_db_download"
    label 'process_low'

    conda "bioconda::eggnog-mapper=2.1.13"
    container "${ workflow.containerEngine == 'singularity' && !task.ext.singularity_pull_docker_container ?
        'https://depot.galaxyproject.org/singularity/eggnog-mapper:2.1.13--pyhdfd78af_2' :
        'biocontainers/eggnog-mapper:2.1.13--pyhdfd78af_2' }"

    output:
    path "eggnog_data", emit: eggnog_data_dir

    when:
    task.ext.when == null || task.ext.when

    script:
    """
    mkdir eggnog_data

    # Download directly from correct URL (download_eggnog_data.py uses wrong host)
    wget -q http://eggnog5.embl.de/download/emapperdb-5.0.2/eggnog.db.gz \\
        -O eggnog_data/eggnog.db.gz
    wget -q http://eggnog5.embl.de/download/emapperdb-5.0.2/eggnog_proteins.dmnd.gz \\
        -O eggnog_data/eggnog_proteins.dmnd.gz
    wget -q http://eggnog5.embl.de/download/emapperdb-5.0.2/eggnog.taxa.tar.gz \\
        -O eggnog_data/eggnog.taxa.tar.gz

    # Decompress
    gunzip -f eggnog_data/eggnog.db.gz
    gunzip -f eggnog_data/eggnog_proteins.dmnd.gz
    tar -xzf eggnog_data/eggnog.taxa.tar.gz -C eggnog_data/ && rm eggnog_data/eggnog.taxa.tar.gz
    """
}
