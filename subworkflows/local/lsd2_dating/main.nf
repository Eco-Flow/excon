//
// --dating_method lsd2: time-calibrate the IQ-TREE species tree with IQ-TREE's LSD2
// instead of ape::chronos, from the same calibration table (--tree_calibrations).
// Three steps, because the IQ-TREE container has no R: check the calibrations and
// write IQ-TREE's date file; date; check the dated tree. See bin/lsd2_dating.R.
//

include { LSD2_PREPARE   } from '../../../modules/local/lsd2_prepare'
include { DATE_TREE_LSD2 } from '../../../modules/local/date_tree_lsd2'
include { LSD2_CHECK     } from '../../../modules/local/lsd2_check'

workflow LSD2_DATING {

    take:
    ch_alignment     // path: the supermatrix the species tree was inferred from
    ch_tree          // path: the rooted species tree (substitutions/site)
    ch_calibrations  // path: calibration table (clade, tips, age_min, age_max)

    main:
    LSD2_PREPARE ( ch_tree, ch_calibrations )

    DATE_TREE_LSD2 (
        ch_alignment,
        ch_tree,
        LSD2_PREPARE.out.dates,
        LSD2_PREPARE.out.outgroup
    )

    LSD2_CHECK (
        DATE_TREE_LSD2.out.nex,
        ch_calibrations,
        ch_tree,
        DATE_TREE_LSD2.out.log
    )

    emit:
    dated_tree   = LSD2_CHECK.out.dated_tree    // path: SpeciesTree_dated.nwk (Myr)
    calibrations = LSD2_CHECK.out.calibrations  // path: dating_calibrations.tsv
    qc           = LSD2_CHECK.out.qc            // path: dating_qc.tsv
    timetree     = DATE_TREE_LSD2.out.nex       // path: IQ-TREE's dated tree (NEXUS)
}
