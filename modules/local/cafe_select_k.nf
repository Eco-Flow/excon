process CAFE_SELECT_K {
    label 'process_single'
    container 'ecoflowucl/chopgo:r-4.3.2_python-3.10_perl-5.38'

    input:
    path result_dirs

    output:
    path "model_selection.tsv", emit: table
    path "best_k.txt",          emit: best_k
    path "k_boundary_warning.txt", emit: warnings, optional: true
    tuple val("${task.process}"), val('python'), val('3.10.13'), emit: versions_python, topic: versions

    script:
    """
    /usr/bin/python3.10 ${projectDir}/bin/cafe_select_k.py ${params.cafe_max_k} ${params.cafe_fixed_lambda != null ? 'fixed_lambda' : ''} ${params.cafe_fixed_alpha != null ? 'fixed_alpha' : ''}
    """

    stub:
    """
    echo -e "k\tn_params\tneg_lnL\tAIC\tBIC\tn_families\tSelected\tNote" > model_selection.tsv
    echo -e "3\t2\t195812.3\t391628.6\tNA\tNA\tBEST\t"                 >> model_selection.tsv
    echo "3" > best_k.txt
    """
}
