#!/usr/bin/env nextflow

log.info """\
=========================================

 EXCON v2.1.0

 -----------------------------------------

 Authors:
   - Chris Wyatt <c.wyatt@ucl.ac.uk>

 -----------------------------------------

 Copyright (c) 2021

 =========================================""".stripIndent()

include { validateParameters; paramsHelp; paramsSummaryLog } from 'plugin/nf-schema'

include { CAFE_PREP } from './modules/local/cafe_prep.nf'
include { CAFE_RUN } from './modules/local/cafe_run.nf'
include { CAFE_MODEL_COMPARE } from './modules/local/cafe_model_compare.nf'
include { CAFE_GO_PREP } from './modules/local/cafe_go_prep.nf'
include { CAFE_GO_RUN } from './modules/local/cafe_go_run.nf'
include { RESCALE_TREE } from './modules/local/rescale_tree.nf'
include { CHROMO_GO } from './modules/local/chromo_go.nf'
include { CAFE_PLOT } from './modules/local/cafe_plot.nf'
include { RENAME_FASTA } from './modules/local/rename_fasta.nf'
include { EGGNOG_DOWNLOAD } from './modules/local/eggnog_download.nf'
include { EGGNOG_TO_GO } from './modules/local/eggnog_to_go.nf'
include { EGGNOG_TO_OG_GO } from './modules/local/eggnog_to_og_go.nf'
include { SUMMARIZE_CHROMO_GO } from './modules/local/sumarize_chromosome_go.nf'

include { NCBIGENOMEDOWNLOAD } from './modules/nf-core/ncbigenomedownload/main.nf'
include { GFFREAD } from './modules/nf-core/gffread/main.nf'
include { BUSCO_BUSCO } from './modules/nf-core/busco/busco/main.nf'
include { AGAT_SPSTATISTICS } from './modules/nf-core/agat/spstatistics/main.nf'
include { AGAT_CONVERTSPGXF2GXF } from './modules/nf-core/agat/convertspgxf2gxf/main.nf'
include { AGAT_SPKEEPLONGESTISOFORM } from './modules/nf-core/agat/spkeeplongestisoform/main.nf'
include { QUAST } from './modules/nf-core/quast/main.nf'
include { GUNZIP } from './modules/nf-core/gunzip/main.nf'
include { ORTHOFINDER as ORTHOFINDER_CAFE } from './modules/nf-core/orthofinder/main.nf'
include { EGGNOGMAPPER } from './modules/nf-core/eggnogmapper/main.nf'

workflow {

   if (params.help) {
      log.info paramsHelp("nextflow run main.nf --input input_file.csv")
      exit 0
   }

   Channel
   .fromPath(params.input)
   .splitCsv()
   .branch {
      ncbi: it.size() == 2
      local: it.size() == 3
   }
   .set { input_type }

   validateParameters()
   log.info paramsSummaryLog(workflow)

   // Write each accession string to its own file, named by sample id
   ch_ncbi = input_type.ncbi
      .collectFile { row -> [ "${row[0]}.txt", row[1] + '\n' ] }
      .map { f -> [ [id: f.baseName], f ] }

   NCBIGENOMEDOWNLOAD (
      ch_ncbi.map { meta, f -> meta },   // val meta
      ch_ncbi.map { meta, f -> f },      // path accessions (now an actual file)
      [],                                // path taxids
      params.groups
   )

   ch_gff = NCBIGENOMEDOWNLOAD.out.gff
      .mix( input_type.local.map { [ [id: it[0]], file(it[2]) ] } )

   ch_fna_raw = NCBIGENOMEDOWNLOAD.out.fna
      .mix( input_type.local.map { [ [id: it[0]], file(it[1]) ] } )

   // Split on .gz, decompress only what needs it
   ch_fna_gz    = ch_fna_raw.filter { meta, fna -> fna.name.endsWith('.gz') }
   ch_fna_plain = ch_fna_raw.filter { meta, fna -> !fna.name.endsWith('.gz') }

   GUNZIP ( ch_fna_gz )
   ch_fna = GUNZIP.out.gunzip.mix( ch_fna_plain )

   // Convert GFF to standard AGAT format then keep longest isoform
   AGAT_CONVERTSPGXF2GXF ( ch_gff )
   AGAT_SPKEEPLONGESTISOFORM ( AGAT_CONVERTSPGXF2GXF.out.output_gff, [] )

   // Join fna + agat gff by meta, then split for GFFREAD's two inputs
   ch_fna_gff = ch_fna.join( AGAT_SPKEEPLONGESTISOFORM.out.gff )

   GFFREAD (
      ch_fna_gff.map { meta, fna, gff -> [ meta, gff ] },
      ch_fna_gff.map { meta, fna, gff -> fna }
   )

   // Remove stop codons from protein fasta
   ch_fasta_for_rename = GFFREAD.out.gffread_fasta.join(
    AGAT_SPKEEPLONGESTISOFORM.out.gff
   )

   RENAME_FASTA (
    ch_fasta_for_rename.map { meta, fasta, gff -> [ meta, fasta ] },
    ch_fasta_for_rename.map { meta, fasta, gff -> [ meta, gff ] }
   )

   // Use renamed fasta for everything downstream
   merge_ch = RENAME_FASTA.out.fasta

   // --- Eggnog GO annotation --- 

   if (params.run_eggnog) {

      if (params.eggnog_data_dir) {
         ch_eggnog_data = channel.fromPath(params.eggnog_data_dir).first()
      } else {
         EGGNOG_DOWNLOAD()
         ch_eggnog_data = EGGNOG_DOWNLOAD.out.eggnog_data_dir.first()
      }

      EGGNOGMAPPER (
         RENAME_FASTA.out.fasta,
         channel.value([ 'diamond', [] ]),
         ch_eggnog_data
      )

      ch_annot_gff = EGGNOGMAPPER.out.annotations.join(
         AGAT_CONVERTSPGXF2GXF.out.output_gff    // all isoforms, not longest only
      )

      EGGNOG_TO_GO (
         ch_annot_gff.map { meta, annot, gff -> [ meta, annot ] },
         ch_annot_gff.map { meta, annot, gff -> [ meta, gff ] }
      )

      ch_go_files = EGGNOG_TO_GO.out.go_file
    	.map { meta, go -> go }
    	.collect()
   }

   // --- Quality stats --- 

   if (params.stats) {
      BUSCO_BUSCO (
         GFFREAD.out.gffread_fasta,
         params.busco_mode,
         params.busco_lineage,
         params.busco_lineages_path ?: [],
         params.busco_config ?: [],
         []
      )

      AGAT_SPSTATISTICS ( AGAT_SPKEEPLONGESTISOFORM.out.gff )

      ch_quast_input = ch_fna
         .join(AGAT_SPKEEPLONGESTISOFORM.out.gff)

      QUAST(
         ch_quast_input.map { meta, fasta, gff -> [ meta, fasta ] },
         ch_quast_input.map { meta, fasta, gff -> [ [], [] ] },   // no reads
         ch_quast_input.map { meta, fasta, gff -> [ meta, gff ] }
      )
   }


    // --- CAFE gene family evolution ---

    if (!params.skip_cafe) {

        ORTHOFINDER_CAFE (
            merge_ch
                .map { meta, fasta -> fasta }
                .collect()
                .map { files -> [ [id: "ortho_cafe"], files ] },
            [[],[]]
        )

        RESCALE_TREE ( ORTHOFINDER_CAFE.out.speciestree )

        CAFE_PREP (
            ORTHOFINDER_CAFE.out.orthologues,
            RESCALE_TREE.out.rescaled_tree
        )

        CAFE_RUN (
            Channel.of( [id: 'gamma'], [id: 'gamma_per_family'] )
                .combine( CAFE_PREP.out.prepared_counts )
                .combine( CAFE_PREP.out.prepared_tree )
                .map { meta, counts, tree -> [ meta, counts, tree ] }
        )

        CAFE_MODEL_COMPARE (
            CAFE_PREP.out.results,
            CAFE_RUN.out.results.filter { meta, res -> meta.id == 'gamma' },
            CAFE_RUN.out.results.filter { meta, res -> meta.id == 'gamma_per_family' }
        )

        // --- Select best model ---
        ch_all_results = CAFE_PREP.out.results
            .map { res -> [ 'base', res ] }
            .mix(
                CAFE_RUN.out.results.map { meta, res -> [ meta.id, res ] }
            )

        ch_best_results = CAFE_MODEL_COMPARE.out.best_model
            .map  { f -> f.text.trim() }
            .combine( ch_all_results )
            .filter { best, model, res -> best == model }
            .map    { best, model, res -> res }

        CAFE_PLOT ( ch_best_results )

        // --- CAFE GO enrichment (requires eggnog) ---

        if (params.run_eggnog) {

            EGGNOG_TO_OG_GO (
                ch_go_files,
                ORTHOFINDER_CAFE.out.orthologues
            )

            CAFE_GO_PREP (
                ch_best_results,
                CAFE_PREP.out.N0_table,
                EGGNOG_TO_OG_GO.out.og_go
            )

            // Flatten pos and neg file lists into individual items
            ch_target_files = CAFE_GO_PREP.out.pos_files
                .mix( CAFE_GO_PREP.out.neg_files )
                .flatten()

            // Flatten background file list into individual items
            ch_bk_files = CAFE_GO_PREP.out.bk_files
                .flatten()

            // Parse manifest - one row per ChopGO job
            ch_manifest = CAFE_GO_PREP.out.manifest
                .splitCsv( sep: '\t', header: false )
                .map { row ->
                    def name = row[0].replaceAll(/\.txt$/, '')
                    tuple( [id: name], row[0], row[1] )
                }

            // Match each manifest row to its target file by filename
            ch_with_target = ch_manifest
                .combine( ch_target_files )
                .filter { meta, target_name, bg_name, file -> file.name == target_name }
                .map    { meta, target_name, bg_name, file -> tuple( meta, file, bg_name ) }

            // Match each row to its background file by filename
            ch_with_bg = ch_with_target
                .combine( ch_bk_files )
                .filter { meta, target_file, bg_name, file -> file.name == bg_name }
                .map    { meta, target_file, bg_name, file -> tuple( meta, target_file, file ) }

            // Add the shared OG_GO file to every job
            ch_go_run_input = ch_with_bg
                .combine( CAFE_GO_PREP.out.og_go.first() )
                .map { meta, target_file, bg_file, og_go ->
                    tuple( meta, target_file, bg_file, og_go )
                }

            CAFE_GO_RUN ( ch_go_run_input )

        } // end if run_eggnog (CAFE GO)

    } // end if !skip_cafe


    // --- Chromosome GO analysis (requires eggnog) ---

    if (params.chromo_go && params.run_eggnog) {

        ch_gff_go = AGAT_SPKEEPLONGESTISOFORM.out.gff
        .join( EGGNOG_TO_GO.out.go_file )
        .map { meta, gff, go -> tuple(meta, gff, go) }

        CHROMO_GO ( ch_gff_go, ORTHOFINDER_CAFE.out.orthologues)

        SUMMARIZE_CHROMO_GO ( CHROMO_GO.out.chromosome_go_filt.mix( CHROMO_GO.out.chromosome_go_unfilt ))

    } // end if chromo_go


   Channel.topic('versions')
    .unique()
    .map { process, tool, version -> "${process}:\n    ${tool}: ${version}" }
    .collectFile(
        name:      'software_versions.yml',
        storeDir:  "${params.outdir}/pipeline_info",
        sort:      true,
        newLine:   true
    )


}

workflow.onComplete {
   println ( workflow.success ? "\nDone!\n" : "Oops... something went wrong" )
}
