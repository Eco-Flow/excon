// Test-only workflow for tests/go_tables/main.nf.test: runs CHROMO_GO (one task per species,
// as the pipeline does), SUMMARIZE_CHROMO_GO and CAFE_GO_RUN on the fixtures in data/, plus an
// independent topGO run per species to check their tables against.

include { CHROMO_GO           } from '../../modules/local/chromo_go.nf'
include { SUMMARIZE_CHROMO_GO } from '../../modules/local/sumarize_chromosome_go.nf'
include { CAFE_GO_RUN         } from '../../modules/local/cafe_go_run.nf'

process INDEPENDENT_TOPGO {
    tag "${meta.id}"
    container 'ecoflowucl/chopgo:r-4.3.2_python-3.10_perl-5.38'

    input:
    path rscript
    tuple val(meta), path(gff), path(go), path(cafe_lists)

    output:
    tuple val(meta), path("independent.tsv"), emit: table

    script:
    """
    Rscript ${rscript} ${gff} ${go} ${cafe_lists}
    """
}

workflow GO_TABLES {
    take:
    rscript
    species      // channel: [ meta, gff, go ]
    orthogroups
    cafe         // channel: [ meta, target, background, og_go ]

    main:
    CHROMO_GO ( species, orthogroups )
    SUMMARIZE_CHROMO_GO ( CHROMO_GO.out.chromosome_go_filt.mix( CHROMO_GO.out.chromosome_go_unfilt ) )
    CAFE_GO_RUN ( cafe )

    // The CAFE GO test set and background go to the independent run of the same species
    cafe_lists = cafe.map { meta, target, background, og_go -> [ meta, [ target, background ] ] }
    INDEPENDENT_TOPGO (
        rscript,
        species.join( cafe_lists, remainder: true ).map { meta, gff, go, lists -> [ meta, gff, go, lists ?: [] ] }
    )

    emit:
    chromo_unfilt = CHROMO_GO.out.chromosome_go_unfilt
    chromo_filt   = CHROMO_GO.out.chromosome_go_filt
    summaries     = SUMMARIZE_CHROMO_GO.out.tables
    cafe          = CAFE_GO_RUN.out.topgo_results
    independent   = INDEPENDENT_TOPGO.out.table
}
