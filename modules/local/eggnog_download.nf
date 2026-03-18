process EGGNOG_DOWNLOAD {
    tag "eggnog_db_download"
    label 'process_low'

    conda "bioconda::eggnog-mapper=2.1.13"
    container 'biocontainers/eggnog-mapper:2.1.13--pyhdfd78af_2'

    output:
    path "eggnog_data", emit: eggnog_data_dir

    when:
    task.ext.when == null || task.ext.when

    script:
    """
    mkdir eggnog_data
    download_eggnog_data.py \\
        --data_dir eggnog_data \\
        -y \\
        -q
    """
}