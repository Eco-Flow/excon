#!/usr/bin/env nextflow

log.info """\
=========================================

 EXCON v2.2.0

 -----------------------------------------

 Authors:
   - Chris Wyatt <c.wyatt@ucl.ac.uk>

 -----------------------------------------

 Copyright (c) 2021

 =========================================""".stripIndent()

include { validateParameters; paramsHelp; paramsSummaryLog } from 'plugin/nf-schema'

include { RESCALE_TREE } from './modules/local/rescale_tree.nf'
include { CAFE_RUN } from './modules/local/cafe_run.nf'
include { CAFE_MODEL_COMPARE } from './modules/local/cafe_model_compare.nf'
include { CAFE_GO_PREP } from './modules/local/cafe_go_prep.nf'
include { CAFE_GO_RUN } from './modules/local/cafe_go_run.nf'
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
include { AGAT_SPKEEPLONGESTISOFORM } from './modules/nf-core/agat/spkeeplongestisoform/main.nf'
include { QUAST } from './modules/nf-core/quast/main.nf'
include { GUNZIP } from './modules/nf-core/gunzip/main.nf'
include { ORTHOFINDER as ORTHOFINDER_CAFE } from './modules/nf-core/orthofinder/main.nf'
include { ORTHOFINDER_V2 as ORTHOFINDER_V2_CAFE } from './modules/local/orthofinder_v2.nf'
include { EGGNOGMAPPER } from './modules/nf-core/eggnogmapper/main.nf'

include { CAFE_PREP } from './modules/local/cafe_prep.nf'
include { CAFE_RUN_K } from './modules/local/cafe_run_k.nf'
include { CAFE_SELECT_K } from './modules/local/cafe_select_k.nf'
include { CAFE_RUN_BEST } from './modules/local/cafe_run_best.nf'
include { CAFE_RUN_LARGE } from './modules/local/cafe_run_large.nf'
include { CAFE_PLOT as CAFE_PLOT_LARGE } from './modules/local/cafe_plot.nf'
include { CAFE_GO_PREP as CAFE_GO_PREP_LARGE } from './modules/local/cafe_go_prep.nf'
include { CAFE_GO_RUN  as CAFE_GO_RUN_LARGE  } from './modules/local/cafe_go_run.nf'

workflow {

   if (params.help) {
      log.info paramsHelp("nextflow run main.nf --input input_file.csv")
      exit 0
   }

   validateParameters()
   log.info paramsSummaryLog(workflow)

   // Whether to skip genome processing (download → AGAT → GFFREAD → RENAME_FASTA → OrthoFinder)
   // Only possible when a pre-computed tree and orthogroups are supplied, AND the user
   // is not requesting EggNOG annotation or genome quality stats (which need the proteins/assemblies).
   def use_precomputed = params.input_tree && params.input_orthogroups
   def needs_genomes   = !use_precomputed || params.run_eggnog || params.stats

   if (needs_genomes && !params.input) {
      error "ERROR: --input (samplesheet CSV) is required when not using pre-computed OrthoFinder results, or when --run_eggnog / --stats is set."
   }

   if (needs_genomes) {

      Channel
      .fromPath(params.input)
      .splitCsv()
      .branch {
         ncbi: it.size() == 2
         local: it.size() == 3
      }
      .set { input_type }

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

      // Keep longest isoform (AGAT sanitises the GFF as part of this step)
      AGAT_SPKEEPLONGESTISOFORM ( ch_gff, [] )

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

   } // end if needs_genomes


   // --- Eggnog GO annotation ---

   if (params.run_eggnog) {

      if (params.eggnog_data_dir) {
         ch_eggnog_data = channel.value(file(params.eggnog_data_dir))
      } else {
         EGGNOG_DOWNLOAD()
         ch_eggnog_data = EGGNOG_DOWNLOAD.out.eggnog_data_dir
      }

      EGGNOGMAPPER (
         RENAME_FASTA.out.fasta,
         channel.value([ 'diamond', [] ]),
         ch_eggnog_data
      )

      ch_annot_gff = EGGNOGMAPPER.out.annotations.join(
         ch_gff    // all isoforms
      )

      EGGNOG_TO_GO (
         ch_annot_gff.map { meta, annot, gff -> [ meta, annot ] },
         ch_annot_gff.map { meta, annot, gff -> [ meta, gff ] }
      )

      ch_go_file_meta = EGGNOG_TO_GO.out.go_file
      ch_go_files = ch_go_file_meta
    	.map { meta, go -> go }
    	.collect()

   } else if (params.predownloaded_gofiles) {

      // User-provided gene-to-GO files (one *.go.txt per species, tab-separated gene_id<TAB>GO:term)
      ch_go_file_meta = Channel.fromPath("${params.predownloaded_gofiles}/*.go.txt")
         .map { file -> [ [id: file.simpleName], file ] }
      ch_go_files = ch_go_file_meta
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

        if (params.input_tree && params.input_orthogroups) {
            ch_speciestree = Channel.fromPath(params.input_tree, checkIfExists: true)
            ch_orthologues = Channel.fromPath(params.input_orthogroups, checkIfExists: true)
        } else if (params.orthofinder_v2) {
            ORTHOFINDER_V2_CAFE (
                merge_ch
                    .map { meta, fasta -> fasta }
                    .collect()
                    .map { files -> [ [id: "ortho_cafe"], files ] }
            )
            ch_speciestree = ORTHOFINDER_V2_CAFE.out.speciestree
            ch_orthologues = ORTHOFINDER_V2_CAFE.out.orthologues
        } else {
            ORTHOFINDER_CAFE (
                merge_ch
                    .map { meta, fasta -> fasta }
                    .collect()
                    .map { files -> [ [id: "ortho_cafe"], files ] },
                [[],[]]
            )
            ch_speciestree = ORTHOFINDER_CAFE.out.speciestree
            ch_orthologues = ORTHOFINDER_CAFE.out.orthologues
        }

        RESCALE_TREE ( ch_speciestree )

        CAFE_PREP (
            ch_orthologues,
            RESCALE_TREE.out.rescaled_tree
        )

        // Run CAFE with fixed lambda on high-differential families filtered out during prep.
        // Only executes when cafe_prep_filtered.R was triggered (attempt > 1) and
        // found families above the differential threshold — otherwise large_counts is empty.
        CAFE_RUN_LARGE (
            CAFE_PREP.out.large_counts,
            CAFE_PREP.out.pruned_tree,
            CAFE_PREP.out.error_model,
            CAFE_PREP.out.lambda.map { f -> f.text.trim() }
        )


        k_values = Channel.of( 1..params.cafe_max_k )

        CAFE_RUN_K (
        CAFE_PREP.out.prepared_counts,   // hog_gene_counts.tsv — possibly filtered
        CAFE_PREP.out.pruned_tree,       // rescaled tree with species names already stripped
        CAFE_PREP.out.error_model,       // Base_error_model.txt — empty file if estimation failed
        k_values                         // each fans out: 1, 2, 3 ... cafe_max_k
        )

        CAFE_SELECT_K(
            CAFE_RUN_K.out.results.map { k, d -> d }.collect()
        )

        // Read the integer out of best_k.txt for passing to CAFE_RUN_BEST
        best_k_ch = CAFE_SELECT_K.out.best_k
           .map { f -> f.text.trim().toInteger() }

        ch_best_uniform = CAFE_RUN_K.out.results
            .combine( best_k_ch )
            .filter { k, dir, best_k -> k == best_k }
            .map    { k, dir, best_k -> dir }

        CAFE_RUN_BEST(
           CAFE_PREP.out.prepared_counts,
           CAFE_PREP.out.pruned_tree,
           CAFE_PREP.out.error_model,
           best_k_ch,
           Channel.of( true )  //Only run poisson here, as we ran without -p earlier
        )

        // Compare uniform vs Poisson at best k, emit the winning directory
        CAFE_MODEL_COMPARE (
            ch_best_uniform,
            CAFE_RUN_BEST.out.results
        )

        ch_best_results = CAFE_MODEL_COMPARE.out.best_model
            .map { f -> f.text.trim() }
            .combine(
                ch_best_uniform.map       { dir -> [ 'uniform', dir ] }
                .mix( CAFE_RUN_BEST.out.results.map { dir -> [ 'poisson', dir ] } )
            )
            .filter { best, model, dir -> best == model }
            .map    { best, model, dir -> dir }


        CAFE_PLOT ( ch_best_results )

        // Plot high-differential families — only runs when CAFE_RUN_LARGE executed
        CAFE_PLOT_LARGE ( CAFE_RUN_LARGE.out.results )


        // --- CAFE GO enrichment (requires eggnog) ---

        if (params.run_eggnog || params.predownloaded_gofiles) {

            EGGNOG_TO_OG_GO (
                ch_go_files,
                ch_orthologues
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
                .combine( CAFE_GO_PREP.out.og_go )
                .map { meta, target_file, bg_file, og_go ->
                    tuple( meta, target_file, bg_file, og_go )
                }

            CAFE_GO_RUN ( ch_go_run_input )


            // --- GO enrichment on high-differential (large) families ---
            // Only fires when CAFE_RUN_LARGE ran (i.e. large_counts was non-empty).
            // Reuses the same EGGNOG_TO_OG_GO output — no extra annotation work needed.

            CAFE_GO_PREP_LARGE (
                CAFE_RUN_LARGE.out.results,
                CAFE_PREP.out.N0_table,
                EGGNOG_TO_OG_GO.out.og_go
            )

            ch_large_target_files = CAFE_GO_PREP_LARGE.out.pos_files
                .mix( CAFE_GO_PREP_LARGE.out.neg_files )
                .flatten()

            ch_large_bk_files = CAFE_GO_PREP_LARGE.out.bk_files
                .flatten()

            ch_large_manifest = CAFE_GO_PREP_LARGE.out.manifest
                .splitCsv( sep: '\t', header: false )
                .map { row ->
                    def name = row[0].replaceAll(/\.txt$/, '')
                    tuple( [id: "large_${name}"], row[0], row[1] )
                }

            ch_large_with_target = ch_large_manifest
                .combine( ch_large_target_files )
                .filter { meta, target_name, bg_name, file -> file.name == target_name }
                .map    { meta, target_name, bg_name, file -> tuple( meta, file, bg_name ) }

            ch_large_with_bg = ch_large_with_target
                .combine( ch_large_bk_files )
                .filter { meta, target_file, bg_name, file -> file.name == bg_name }
                .map    { meta, target_file, bg_name, file -> tuple( meta, target_file, file ) }

            ch_large_go_run_input = ch_large_with_bg
                .combine( CAFE_GO_PREP_LARGE.out.og_go )
                .map { meta, target_file, bg_file, og_go ->
                    tuple( meta, target_file, bg_file, og_go )
                }

            CAFE_GO_RUN_LARGE ( ch_large_go_run_input )

        } // end if run_eggnog / predownloaded_gofiles (CAFE GO)

    } // end if !skip_cafe


    // --- Chromosome GO analysis (requires eggnog) ---

    if (params.chromo_go && (params.run_eggnog || params.predownloaded_gofiles)) {

        ch_gff_go = AGAT_SPKEEPLONGESTISOFORM.out.gff
        .join( ch_go_file_meta )
        .map { meta, gff, go -> tuple(meta, gff, go) }

        CHROMO_GO ( ch_gff_go, ch_orthologues)

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
