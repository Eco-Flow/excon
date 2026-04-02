process CAFE_SELECT_K {
    label 'process_single'
    container 'ecoflowucl/chopgo:r-4.3.2_python-3.10_perl-5.38'

    input:
    path result_dirs

    output:
    path "model_selection.tsv", emit: table
    path "best_k.txt",          emit: best_k
    tuple val("${task.process}"), val('python'), val('3.10.13'), emit: versions_python, topic: versions

    script:
    """
    /usr/bin/python3.10 ${projectDir}/bin/cafe_select_k.py
    """

    stub:
    """
    echo -e "k\tn_params\tneg_lnL\tAIC\tBIC\tn_families\tSelected" > model_selection.tsv
    echo -e "3\t2\t195812.3\t391628.6\tNA\tNA\tBEST"               >> model_selection.tsv
    echo "3" > best_k.txt
    """
}
