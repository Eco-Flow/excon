// Test-only workflow for tests/go_tables/main.nf.test: runs CHROMO_GO and CAFE_GO_RUN on
// the BRAKER-style fixture in data/, plus an independent topGO run to check their tables against.

include { CHROMO_GO   } from '../../modules/local/chromo_go.nf'
include { CAFE_GO_RUN } from '../../modules/local/cafe_go_run.nf'

process INDEPENDENT_TOPGO {
    container 'ecoflowucl/chopgo:r-4.3.2_python-3.10_perl-5.38'

    input:
    path rscript
    path gff
    path go
    path target
    path background

    output:
    path "independent.tsv", emit: table

    script:
    """
    Rscript ${rscript} ${gff} ${go} ${target} ${background}
    """
}

workflow GO_TABLES {
    take:
    rscript
    gff
    go
    orthogroups
    target
    background
    og_go

    main:
    CHROMO_GO ( gff.combine(go).map { g, t -> [ [id: 'Testsp'], g, t ] }, orthogroups )
    CAFE_GO_RUN ( target.combine(background).combine(og_go).map { t, b, o -> [ [id: 'Testsp'], t, b, o ] } )
    INDEPENDENT_TOPGO ( rscript, gff, go, target, background )

    emit:
    chromo_unfilt = CHROMO_GO.out.chromosome_go_unfilt
    chromo_filt   = CHROMO_GO.out.chromosome_go_filt
    cafe          = CAFE_GO_RUN.out.topgo_results
    independent   = INDEPENDENT_TOPGO.out.table
}
