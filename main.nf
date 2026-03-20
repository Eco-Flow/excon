#!/usr/bin/env nextflow

log.info """\
 =========================================

 EXCON (v2.0.0)

 -----------------------------------------

 Authors:
   - Chris Wyatt <c.wyatt@ucl.ac.uk>

 -----------------------------------------

 Copyright (c) 2021

 =========================================""".stripIndent()

include { validateParameters; paramsHelp; paramsSummaryLog } from 'plugin/nf-validation'

include { CAFE } from './modules/local/cafe_with_retry.nf'
include { RESCALE_TREE } from './modules/local/rescale_tree.nf'
include { CHROMO_GO } from './modules/local/chromo_go.nf'
include { CAFE_GO } from './modules/local/cafe_go.nf'
include { CAFE_PLOT } from './modules/local/cafe_plot.nf'
include { FILTER_FASTA } from './modules/local/filter_fasta'
include { EGGNOG_DOWNLOAD } from './modules/local/eggnog_download.nf'
include { EGGNOG_TO_GO } from './modules/local/eggnog_to_go.nf'
include { EGGNOG_TO_OG_GO } from './modules/local/eggnog_to_og_go.nf'

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
include { CUSTOM_DUMPSOFTWAREVERSIONS } from './modules/nf-core/custom/dumpsoftwareversions/main.nf'

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

   ch_versions = Channel.empty()

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
   FILTER_FASTA ( GFFREAD.out.gffread_fasta )
   merge_ch = FILTER_FASTA.out

   // --- Eggnog GO annotation --- 

   if (params.run_eggnog) {

      if (params.eggnog_data_dir) {
         ch_eggnog_data = channel.fromPath(params.eggnog_data_dir).first()
      } else {
         EGGNOG_DOWNLOAD()
         ch_eggnog_data = EGGNOG_DOWNLOAD.out.eggnog_data_dir.first()
      }

      EGGNOGMAPPER (
         FILTER_FASTA.out,
         channel.value([ 'diamond', [] ]),
         ch_eggnog_data
      )
      ch_versions = ch_versions.mix(EGGNOGMAPPER.out.versions_eggnogmapper.first())

      EGGNOG_TO_GO ( EGGNOGMAPPER.out.annotations )

      ch_go_files = EGGNOG_TO_GO.out.go_file
         .map { meta, go -> go }
         .collect()
   }

   // --- Quality stats --- 

   if (params.stats) {
      BUSCO_BUSCO (
         GFFREAD.out.proteins_busco,
         params.busco_mode,
         params.busco_lineage,
         params.busco_lineages_path ?: [],
         params.busco_config ?: []
      )
      ch_versions = ch_versions.mix(BUSCO_BUSCO.out.versions.first())

      AGAT_SPSTATISTICS ( GFFREAD.out.gffs_agat )
      ch_versions = ch_versions.mix(AGAT_SPSTATISTICS.out.versions.first())

      QUAST ( GFFREAD.out.fasta_quast, GFFREAD.out.gffs_agat )
      ch_versions = ch_versions.mix(QUAST.out.versions.first())
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

      CAFE ( ORTHOFINDER_CAFE.out.orthologues, RESCALE_TREE.out.rescaled_tree )
      ch_versions = ch_versions.mix(CAFE.out.versions)

      CAFE_PLOT ( CAFE.out.result )
      ch_versions = ch_versions.mix(CAFE_PLOT.out.versions)

      // -- CAFE GO enrichment (requires eggnog) ---

      if (params.run_eggnog) {
         EGGNOG_TO_OG_GO (
            ch_go_files,
            ORTHOFINDER_CAFE.out.orthologues
         )

         CAFE_GO (
            CAFE.out.result,
            CAFE.out.N0_table,
            EGGNOG_TO_OG_GO.out.og_go
         )
         ch_versions = ch_versions.mix(CAFE_GO.out.versions.first())
      }

      // --- Chromosome GO analysis (requires eggnog) ---

      if (params.chromo_go && params.run_eggnog) {
         CHROMO_GO (
            AGAT_SPKEEPLONGESTISOFORM.out.gff.map { meta, gff -> gff }.collect(),
            ch_go_files,
            ORTHOFINDER_CAFE.out.orthologues
         )
         ch_versions = ch_versions.mix(CHROMO_GO.out.versions)
      }
   }

   CUSTOM_DUMPSOFTWAREVERSIONS ( ch_versions.collectFile(name: 'collated_versions.yml') )
}

workflow.onComplete {
   println ( workflow.success ? "\nDone!\n" : "Oops... something went wrong" )
}
