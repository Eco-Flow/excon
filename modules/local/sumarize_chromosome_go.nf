process SUMMARIZE_CHROMO_GO {
    tag "summarize_go"

    label 'process_medium'
    label 'process_single'

    container 'ecoflowucl/chopgo:r-4.3.2_python-3.10_perl-5.38'

    input:
    path res_tabs

    output:
    path "*.pdf", emit: plots
    path "*.csv", emit: tables
    path "*.txt", emit: reports

    script:
    """
    # Copy all input files into working dir (think we can remove)
    # cp ${res_tabs} .

    summarize_chromosome_go.R
    """
}
