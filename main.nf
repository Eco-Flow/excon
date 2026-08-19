#!/usr/bin/env nextflow

include { validateParameters; paramsHelp; paramsSummaryLog } from 'plugin/nf-schema'

include { RESCALE_TREE } from './modules/local/rescale_tree.nf'
include { DATE_TREE } from './modules/local/date_tree.nf'
include { PRUNE_TREE } from './modules/local/prune_tree.nf'
include { CAFE_RUN } from './modules/local/cafe_run.nf'
include { CAFE_MODEL_COMPARE } from './modules/local/cafe_model_compare.nf'
include { CAFE_GO_PREP } from './modules/local/cafe_go_prep.nf'
include { CAFE_GO_RUN } from './modules/local/cafe_go_run.nf'
include { CHROMO_GO } from './modules/local/chromo_go.nf'
include { CAFE_PLOT } from './modules/local/cafe_plot.nf'
include { CAFE_NODE_GUIDE } from './modules/local/cafe_node_guide.nf'
include { CAFE_SIG_FAMILIES } from './modules/local/cafe_sig_families.nf'
include { CAFE_PLOT_ALTVIZ } from './modules/local/cafe_plot_altviz.nf'
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
include { ORTHOFINDER_BLAST as ORTHOFINDER_BLAST_CAFE } from './modules/local/orthofinder_blast.nf'
include { ORTHOFINDER_PHYLO as ORTHOFINDER_PHYLO_CAFE } from './modules/local/orthofinder_phylo.nf'
include { ORTHOFINDER_V2 as ORTHOFINDER_V2_CAFE } from './modules/local/orthofinder_v2.nf'
include { EXTRACT_SINGLE_COPY } from './modules/local/extract_single_copy.nf'
include { ALIGN_SINGLE_COPY } from './modules/local/align_single_copy.nf'
include { CONCAT_SINGLE_COPY } from './modules/local/concat_single_copy.nf'
include { ROOT_TREE } from './modules/local/root_tree.nf'
include { IQTREE as IQTREE_SPECIES_TREE } from './modules/nf-core/iqtree/main.nf'
include { EGGNOGMAPPER } from './modules/nf-core/eggnogmapper/main.nf'

include { CAFE_PREP } from './modules/local/cafe_prep.nf'
include { CAFE_RUN_K } from './modules/local/cafe_run_k.nf'
include { CAFE_SELECT_K } from './modules/local/cafe_select_k.nf'
include { CAFE_RUN_BEST } from './modules/local/cafe_run_best.nf'
include { CAFE_RUN_LARGE } from './modules/local/cafe_run_large.nf'
include { MERGE_CAFE_LARGE_RESULTS } from './modules/local/merge_cafe_large_results.nf'
include { CAFE_PLOT as CAFE_PLOT_LARGE } from './modules/local/cafe_plot.nf'
include { CAFE_SIG_FAMILIES as CAFE_SIG_FAMILIES_LARGE } from './modules/local/cafe_sig_families.nf'
include { COMBINE_CAFE_SIG_FAMILIES } from './modules/local/combine_cafe_sig_families.nf'
include { CAFE_GO_PREP as CAFE_GO_PREP_LARGE } from './modules/local/cafe_go_prep.nf'
include { CAFE_GO_RUN  as CAFE_GO_RUN_LARGE  } from './modules/local/cafe_go_run.nf'
include { SUMMARIZE_CAFE_GO }                                                  from './modules/local/summarize_cafe_go.nf'
include { SUMMARIZE_CAFE_GO as SUMMARIZE_CAFE_GO_LARGE }                       from './modules/local/summarize_cafe_go.nf'
include { PLOT_CAFE_GO }                                                        from './modules/local/plot_cafe_go.nf'
include { PLOT_CAFE_GO      as PLOT_CAFE_GO_LARGE }                            from './modules/local/plot_cafe_go.nf'
include { OG_ANNOTATION_SUMMARY } from './modules/local/og_annotation_summary.nf'

// Parameter types. Nextflow 26 passes every --flag on the command line as a
// String unless the parameter is declared with a type here; defaults stay in
// nextflow.config.
params {
    outdir                    : String
    input                     : String
    chromo_go                 : Boolean
    chromo_go_max_chroms      : Integer
    go_cutoff                 : Float
    go_type                   : String
    go_max_plot               : Integer
    go_algo                   : String
    forks                     : Integer
    clean                     : Boolean
    config_profile_description : String
    config_profile_contact     : String
    config_profile_url         : String
    custom_config             : String
    max_memory                : String
    max_cpus                  : Integer
    max_time                  : String
    trace_report_suffix       : String
    help                      : Boolean
    publish_dir_mode          : String
    groups                    : String
    ncbi_max_forks            : Integer
    stats                     : Boolean
    busco_mode                : String
    busco_lineage             : String
    busco_lineages_path       : String
    busco_config              : String
    internal_stop_action      : String
    input_tree                : String
    input_orthogroups         : String
    orthofinder_blast_results : String
    orthofinder_results       : String
    proteome_dir              : String
    orthofinder_v2            : Boolean
    orthofinder_method        : String
    orthofinder_msa_prog      : String
    orthofinder_search        : String
    orthofinder_tree          : String
    iqtree_species_tree       : Boolean
    iqtree_outgroup           : String
    iqtree_args               : String
    iqtree_partition_model    : String
    iqtree_partition_file     : String
    cafe_clade                : String
    cafe_species              : String
    tree_calibrations         : String
    chronos_model             : String
    chronos_lambda            : Float
    chronos_rate_categories   : Integer
    tree_scale_factor         : Integer
    input_tree_is_dated       : Boolean
    cafe_zero_root            : Boolean
    skip_cafe                 : Boolean
    cafe_max_differential     : Integer
    cafe_filter_first         : Boolean
    cafe_max_k                : Integer
    orthofinder_msa_dir       : String
    orthofinder_genetree_dir  : String
    cafe_focus_clades         : String
    run_eggnog                : Boolean
    eggnog_data_dir           : String
    eggnog_target_taxa        : String
    eggnog_tax_scope          : String
    eggnog_evalue             : Float
    eggnog_score              : Float
    eggnog_pident             : Float
    eggnog_query_cover        : Float
    eggnog_subject_cover      : Float
    eggnog_rep_species        : String
    predownloaded_gofiles     : String
}

workflow {

   log.info """\
   =========================================
   
    EXCON v${workflow.manifest.version}
   
    -----------------------------------------
   
    Authors:
      - Chris Wyatt <c.wyatt@ucl.ac.uk>
   
    -----------------------------------------
   
    Copyright (c) 2021
   
    =========================================""".stripIndent()

   if (params.help) {
      log.info paramsHelp(command: "nextflow run main.nf --input input_file.csv")
      exit 0
   }

   validateParameters()
   log.info paramsSummaryLog(workflow)

   // Whether to skip genome processing (download → AGAT → GFFREAD → RENAME_FASTA → OrthoFinder)
   // Only possible when a pre-computed tree and orthogroups are supplied, AND the user
   // is not requesting EggNOG annotation or genome quality stats (which need the proteins/assemblies).
   def use_precomputed = (params.input_tree && params.input_orthogroups) || params.orthofinder_blast_results || params.orthofinder_results
   def needs_genomes   = !use_precomputed || params.run_eggnog || params.stats

   // --input_tree and --input_orthogroups only skip OrthoFinder together: the tree
   // says which species to analyse, the table supplies the gene counts.
   if (params.input_tree && !params.input_orthogroups && !params.orthofinder_results) {
      error "ERROR: --input_tree also needs --input_orthogroups (the Orthogroups.tsv or N0.tsv holding the gene counts). Supplied alone it cannot skip OrthoFinder, so --input is required as well."
   }
   if (params.input_orthogroups && !params.input_tree && !params.orthofinder_results) {
      error "ERROR: --input_orthogroups also needs --input_tree (the species tree for those gene counts)."
   }

   if (needs_genomes && !params.input) {
      error "ERROR: --input (samplesheet CSV) is required when not using pre-computed OrthoFinder results, or when --run_eggnog / --stats is set."
   }

   if (params.cafe_clade && params.cafe_species) {
      error "ERROR: give either --cafe_clade or --cafe_species, not both."
   }

   if (!(params.internal_stop_action in ['strip', 'drop', 'longest_orf'])) {
      error "ERROR: --internal_stop_action must be 'strip', 'drop' or 'longest_orf', got '${params.internal_stop_action}'."
   }

   if (params.proteome_dir && !params.orthofinder_results) {
      error "ERROR: --proteome_dir is only used alongside --orthofinder_results."
   }

   // The supermatrix is built from OrthoFinder's single-copy orthogroups, so there is
   // nothing to build it from when a pre-computed tree replaces the OrthoFinder run.
   if (params.iqtree_species_tree) {
      if (params.input_tree && params.input_orthogroups && !params.orthofinder_results) {
         error "ERROR: --iqtree_species_tree needs OrthoFinder output. Either let OrthoFinder run, or point --orthofinder_results at a completed run; --input_tree/--input_orthogroups alone skip OrthoFinder entirely."
      }
      if (params.orthofinder_results && !params.proteome_dir) {
         error "ERROR: --iqtree_species_tree with --orthofinder_results also needs --proteome_dir, holding the proteomes whose gene IDs appear in Orthogroups.tsv (the *.clean.fasta files written by RENAME_FASTA)."
      }
      // '--iqtree_outgroup null' on the command line is the literal string "null",
      // not an unset value, and would only fail once ROOT_TREE runs.
      if (params.iqtree_outgroup?.toLowerCase() in ['null', 'none', 'false']) {
         error "ERROR: --iqtree_outgroup was given the literal value '${params.iqtree_outgroup}'. Omit the option entirely to midpoint-root the tree."
      }
      if (params.iqtree_args?.toLowerCase() in ['null', 'none', 'false']) {
         error "ERROR: --iqtree_args was given the literal value '${params.iqtree_args}', which would be passed to IQ-TREE verbatim. Omit the option entirely."
      }
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

      // OG functional annotation summary — one row per OG with representative gene description
      ch_annot_files = EGGNOGMAPPER.out.annotations
          .map { meta, annot -> annot }
          .collect()

   } else if (params.predownloaded_gofiles) {

      // User-provided gene-to-GO files (one *.go.txt per species, tab-separated gene_id<TAB>GO:term)
      ch_go_file_meta = Channel.fromPath("${params.predownloaded_gofiles}/*.go.txt")
         .map { file -> [ [id: file.simpleName], file ] }
      ch_go_files = ch_go_file_meta
         .map { meta, go -> go }
         .collect()

      ch_annot_files = Channel.empty()

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

    // OrthoFinder outputs are needed for CAFE and/or for the IQ-TREE species tree,
    // so this stage runs whenever either is requested.
    if (!params.skip_cafe || params.iqtree_species_tree) {

        ch_orthofinder_dir = Channel.empty()

        if (params.orthofinder_results) {
            // Reuse a completed OrthoFinder run rather than repeating it.
            ch_orthofinder_dir = Channel.fromPath(params.orthofinder_results, type: 'dir', checkIfExists: true)
                .map { d -> [ [id: "ortho_cafe"], d ] }
            ch_speciestree = Channel.fromPath("${params.orthofinder_results}/Species_Tree/SpeciesTree_rooted_node_labels.txt", checkIfExists: true)
            ch_orthologues = params.input_orthogroups ?
                Channel.fromPath(params.input_orthogroups, checkIfExists: true) :
                Channel.fromPath("${params.orthofinder_results}/Phylogenetic_Hierarchical_Orthogroups/N0.tsv", checkIfExists: true)
        } else if (params.input_tree && params.input_orthogroups) {
            ch_speciestree = Channel.fromPath(params.input_tree, checkIfExists: true)
            ch_orthologues = Channel.fromPath(params.input_orthogroups, checkIfExists: true)
        } else if (params.orthofinder_v2) {
            ORTHOFINDER_V2_CAFE (
                merge_ch
                    .map { meta, fasta -> fasta }
                    .collect()
                    .map { files -> [ [id: "ortho_cafe"], files ] }
            )
            ch_speciestree     = ORTHOFINDER_V2_CAFE.out.speciestree
            ch_orthologues     = ORTHOFINDER_V2_CAFE.out.orthologues
            ch_orthofinder_dir = ORTHOFINDER_V2_CAFE.out.orthofinder
        } else {
            // Stage 1: reciprocal DIAMOND blast (CPU-heavy, parallelisable)
            // Stage 2: orthogroup inference + phylogeny (different resource profile)
            if (params.orthofinder_blast_results) {
                // Resume from a pre-computed blast WorkingDirectory (e.g. from a previous run)
                ch_blast_wd = Channel.fromPath(params.orthofinder_blast_results, checkIfExists: true)
                    .map { d -> [ [id: "ortho_cafe"], d ] }
            } else {
                ORTHOFINDER_BLAST_CAFE (
                    merge_ch
                        .map { meta, fasta -> fasta }
                        .collect()
                        .map { files -> [ [id: "ortho_cafe"], files ] }
                )
                ch_blast_wd = ORTHOFINDER_BLAST_CAFE.out.working_dir
            }

            ORTHOFINDER_PHYLO_CAFE ( ch_blast_wd )
            ch_speciestree     = ORTHOFINDER_PHYLO_CAFE.out.speciestree
            ch_orthologues     = ORTHOFINDER_PHYLO_CAFE.out.orthologues
            ch_orthofinder_dir = ORTHOFINDER_PHYLO_CAFE.out.orthofinder
        }

        // Optionally replace the OrthoFinder species tree with an IQ-TREE2 tree
        // inferred from a concatenated alignment of the single-copy orthogroups.
        // The sequences are pulled from the per-species proteomes and aligned here
        // rather than reusing OrthoFinder's own alignments, so this works whatever
        // --orthofinder_method is set to.
        if (params.iqtree_species_tree) {
            ch_proteomes = params.proteome_dir ?
                Channel.fromPath("${params.proteome_dir}/*", checkIfExists: true).collect() :
                merge_ch.map { meta, fasta -> fasta }.collect()

            EXTRACT_SINGLE_COPY ( ch_orthofinder_dir.combine(ch_proteomes.map { files -> [ files ] }) )

            ALIGN_SINGLE_COPY ( EXTRACT_SINGLE_COPY.out.orthogroups )

            CONCAT_SINGLE_COPY ( ch_orthofinder_dir.join(ALIGN_SINGLE_COPY.out.alignments) )

            // A partition file from a previous run carries that run's per-partition
            // models, so ModelFinder does not have to be repeated.
            ch_partitions = params.iqtree_partition_file ?
                Channel.fromPath(params.iqtree_partition_file, checkIfExists: true) :
                CONCAT_SINGLE_COPY.out.partitions.map { meta, parts -> parts }

            IQTREE_SPECIES_TREE (
                CONCAT_SINGLE_COPY.out.alignment.map { meta, aln -> [ meta, aln, [] ] },
                [], [], [], [],
                ch_partitions,
                [], [], [], [], [], [], []
            )

            ROOT_TREE ( IQTREE_SPECIES_TREE.out.phylogeny )
            ch_speciestree = ROOT_TREE.out.tree
        }
    }

    if (!params.skip_cafe) {

        // A dated, time-calibrated tree is passed to CAFE unchanged (no rescaling);
        // an OrthoFinder substitution tree is scaled first to avoid CAFE5 precision issues.
        if (params.input_tree_is_dated) {
            ch_tree_for_prep = ch_speciestree
        } else if (params.tree_calibrations) {
            // Calibrate here instead, so the tree reaching CAFE5 is on a time axis
            // and lambda is per Myr. Works for the OrthoFinder tree and the
            // IQ-TREE one alike, since both arrive on ch_speciestree.
            DATE_TREE (
                ch_speciestree,
                Channel.fromPath(params.tree_calibrations, checkIfExists: true)
            )
            ch_tree_for_prep = DATE_TREE.out.dated_tree
        } else {
            RESCALE_TREE ( ch_speciestree )
            ch_tree_for_prep = RESCALE_TREE.out.rescaled_tree
        }

        // Restricting to a clade happens last, so dating uses the full tree and
        // calibrations may reference species outside the subset. Pruning preserves
        // node ages, so a calibrated tree stays calibrated.
        if (params.cafe_clade || params.cafe_species) {
            // --cafe_species may be a path to a file of tip names rather than an inline
            // comma-separated list. Under a container profile that file only exists on
            // the launch host, so it must be staged as a proper input or prune_tree.R's
            // file.exists() check fails inside the task's work directory.
            ch_cafe_species_file = (params.cafe_species && file(params.cafe_species).exists()) ?
                Channel.fromPath(params.cafe_species, checkIfExists: true) :
                Channel.fromPath("${projectDir}/assets/NO_FILE")

            PRUNE_TREE ( ch_tree_for_prep, ch_cafe_species_file )
            ch_tree_for_prep = PRUNE_TREE.out.tree
        }

        CAFE_PREP (
            ch_orthologues,
            ch_tree_for_prep
        )

        // High-differential families filtered out during prep are run one at a time,
        // each fitting its own lambda, rather than forcing one lambda to explain every
        // excluded family at once (which never converged). Split with splitCsv/
        // collectFile — entirely inside Nextflow's own dataflow engine — rather than
        // an external process writing one file per family and having Nextflow glob
        // that directory afterwards: on a networked work directory (GPFS etc.) that
        // glob can race the writes and see only some of the files. splitCsv/
        // collectFile emit each family as its own channel item directly, so there is
        // no separate directory listing step to race. Named by HOG id so CAFE_RUN_LARGE
        // task tags and output directories stay readable. Only produces items when
        // cafe_prep_filtered.R was triggered (attempt > 1) and found families above the
        // differential threshold — otherwise large_counts is empty and nothing
        // downstream of it fires.
        ch_large_family_tables = CAFE_PREP.out.large_counts
            .splitCsv( header: true, sep: '\t' )
            .collectFile { row ->
                def hog    = (row.HOG as String).replaceAll(/[^A-Za-z0-9_.-]/, '_')
                def cols   = row.keySet() as List
                def header = cols.join('\t')
                def line   = cols.collect { row[it] }.join('\t')
                [ "${hog}.tsv", "${header}\n${line}\n" ]
            }

        // Explicit combine() rather than passing cafe_tree/error_model as separate
        // positional channels — makes their reuse across every hog_counts item
        // unambiguous rather than relying on Nextflow's implicit broadcast pairing.
        CAFE_RUN_LARGE (
            ch_large_family_tables
                .combine( CAFE_PREP.out.cafe_tree )
                .combine( CAFE_PREP.out.error_model )
        )

        // Stitch the many single-family runs back into one CAFE5-shaped directory so
        // the existing downstream consumers (cafeplotter, CAFE_GO_PREP_LARGE, CAFE_SIG_
        // FAMILIES_LARGE below) can read it exactly like an ordinary CAFE5 result.
        // collect() emits a single empty list even when the source channel had zero
        // items (e.g. no large families at all, or none converged) — filtered out so
        // MERGE_CAFE_LARGE_RESULTS only runs when there is actually something to merge.
        MERGE_CAFE_LARGE_RESULTS (
            CAFE_RUN_LARGE.out.results.collect().filter { it.size() > 0 }
        )


        k_values = Channel.of( 1..params.cafe_max_k )

        CAFE_RUN_K (
        CAFE_PREP.out.prepared_counts,   // hog_gene_counts.tsv — possibly filtered
        CAFE_PREP.out.cafe_tree,         // ultrametric CAFE tree (dated & unchanged, or chronoMPL-scaled)
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
           CAFE_PREP.out.cafe_tree,
           CAFE_PREP.out.error_model,
           best_k_ch,
           Channel.of( true )  //Only run poisson here, as we ran without -p earlier
        )

        // Compare uniform vs Poisson at best k, emit the winning directory
        CAFE_MODEL_COMPARE (
            ch_best_uniform,
            CAFE_RUN_BEST.out.results
        )

        ch_best_results = CAFE_MODEL_COMPARE.out.best_results

        CAFE_PLOT ( ch_best_results )
        CAFE_NODE_GUIDE ( ch_best_results )

        // Per-node significant-family tables, plus (optionally) the alignments and
        // gene trees of families significantly expanded at the focus node(s).
        ch_msa_dir = params.orthofinder_msa_dir ?
            Channel.fromPath(params.orthofinder_msa_dir, type: 'dir', checkIfExists: true) :
            Channel.fromPath("${projectDir}/assets/NO_FILE")
        ch_genetree_dir = params.orthofinder_genetree_dir ?
            Channel.fromPath(params.orthofinder_genetree_dir, type: 'dir', checkIfExists: true) :
            Channel.fromPath("${projectDir}/assets/NO_FILE")

        CAFE_SIG_FAMILIES (
            ch_best_results,
            CAFE_PREP.out.N0_table,
            ch_msa_dir,
            ch_genetree_dir
        )

        // Plot high-differential families — only runs when at least one of them
        // converged on its own lambda (converged.txt is an optional output; if
        // absent the channel is empty and downstream steps are silently skipped)
        ch_large_results_ok = MERGE_CAFE_LARGE_RESULTS.out.converged
            .combine( MERGE_CAFE_LARGE_RESULTS.out.results )
            .map { flag, dir -> dir }

        CAFE_PLOT_LARGE ( ch_large_results_ok )

        // Same significant-family extraction as the main model, run on the merged
        // large-family results — falls back to the NO_FILE placeholder (which
        // cafe_sig_families.R treats as "no CAFE tables found", writing empty
        // tables) when no families were large enough to need this track at all,
        // so COMBINE_CAFE_SIG_FAMILIES below always has something to read.
        ch_large_results_for_sig = ch_large_results_ok
            .ifEmpty( file("${projectDir}/assets/NO_FILE") )

        CAFE_SIG_FAMILIES_LARGE (
            ch_large_results_for_sig,
            CAFE_PREP.out.N0_table,
            ch_msa_dir,
            ch_genetree_dir
        )

        // One combined report spanning every orthogroup CAFE5 could fit at all —
        // the shared k/lambda model's families plus the large-differential ones
        // fit individually — tagged by which track produced each row. Node ids
        // line up between the two: both stem from the same CAFE_PREP.out.cafe_tree,
        // and CAFE5 numbers internal nodes from the tree topology alone.
        COMBINE_CAFE_SIG_FAMILIES (
            CAFE_SIG_FAMILIES.out.changes,
            CAFE_SIG_FAMILIES.out.sig_changes,
            CAFE_SIG_FAMILIES_LARGE.out.changes,
            CAFE_SIG_FAMILIES_LARGE.out.sig_changes
        )


        if (params.run_eggnog) {
        OG_ANNOTATION_SUMMARY (
          ch_annot_files,
          ch_orthologues,
          params.eggnog_rep_species ?: ""
        )
        }

        // --- CAFE GO enrichment ---

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

            CAFE_PLOT_ALTVIZ (
                ch_best_results,
                CAFE_GO_PREP.out.cafe_summary
            )

            CAFE_GO_RUN ( ch_go_run_input )

            SUMMARIZE_CAFE_GO (
                CAFE_GO_RUN.out.topgo_results
                    .map { meta, f -> f }
                    .collect()
                    .map { files -> tuple( "cafe_go", files ) },
                CAFE_GO_RUN.out.pdfs.map { meta, f -> f }.flatten().collect(),
                CAFE_GO_RUN.out.svgs.map { meta, f -> f }.flatten().collect()
            )

            PLOT_CAFE_GO (
                SUMMARIZE_CAFE_GO.out.pos_tsv
                    .join( SUMMARIZE_CAFE_GO.out.neg_tsv )
            )

            // --- GO enrichment on high-differential (large) families ---
            // Only fires when CAFE_RUN_LARGE ran (i.e. large_counts was non-empty).
            // Reuses the same EGGNOG_TO_OG_GO output — no extra annotation work needed.

            CAFE_GO_PREP_LARGE (
                ch_large_results_ok,
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

            SUMMARIZE_CAFE_GO_LARGE (
                CAFE_GO_RUN_LARGE.out.topgo_results
                    .map { meta, f -> f }
                    .collect()
                    .map { files -> tuple( "cafe_go_large", files ) },
                CAFE_GO_RUN_LARGE.out.pdfs.map { meta, f -> f }.flatten().collect(),
                CAFE_GO_RUN_LARGE.out.svgs.map { meta, f -> f }.flatten().collect()
            )

            PLOT_CAFE_GO_LARGE (
                SUMMARIZE_CAFE_GO_LARGE.out.pos_tsv
                    .join( SUMMARIZE_CAFE_GO_LARGE.out.neg_tsv )
            )

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


   workflow.onComplete {
      completionSummary()
   }

}

def completionSummary() {
   println ( workflow.success ? "\nDone!\n" : "Oops... something went wrong" )
}
