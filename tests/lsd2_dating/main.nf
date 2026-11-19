// Test-only workflow for tests/lsd2_dating/main.nf.test: the three LSD2 dating steps,
// chained as main.nf chains them for --dating_method lsd2.

include { LSD2_PREPARE; DATE_TREE_LSD2; LSD2_CHECK } from '../../modules/local/date_tree_lsd2.nf'

workflow LSD2_DATING {
    take:
    alignment
    tree
    calibrations

    main:
    LSD2_PREPARE ( tree, calibrations )
    DATE_TREE_LSD2 ( alignment, tree, LSD2_PREPARE.out.dates, LSD2_PREPARE.out.outgroup )
    LSD2_CHECK ( DATE_TREE_LSD2.out.nex, calibrations, tree, DATE_TREE_LSD2.out.log )

    emit:
    dated_tree   = LSD2_CHECK.out.dated_tree
    calibrations = LSD2_CHECK.out.calibrations
}
