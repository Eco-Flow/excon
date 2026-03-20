#!/usr/bin/env nextflow

log.info """\
 =========================================

 EXCON (v3.0)

 -----------------------------------------

 Authors:
   - Chris Wyatt <c.wyatt@ucl.ac.uk>
   - Simon Murray

 -----------------------------------------

 Copyright (c) 2021

 =========================================""".stripIndent()

include { validateParameters; paramsHelp; paramsSummaryLog } from 'plugin/nf-validation'

include { GET_DATA } from './modules/local/getdata.nf'
include { GO_ASSIGN } from './modules/local/go_assign.nf'
include { GO_EXPANSION  } from './modules/local/go_expansion.nf'
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
include { ORTHOFINDER as ORTHOFINDER_GO} from './modules/nf-core/orthofinder/main.nf'
include { ORTHOFINDER as ORTHOFINDER_CAFE } from './modules/nf-core/orthofinder/main.nf'
include { EGGNOGMAPPER    } from './modules/nf-core/eggnogmapper/main.nf'
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

   //Make a channel for version outputs:
   ch_versions = Channel.empty()

   // Validate input parameters
   validateParameters()

   // Print summary of supplied parameters
   log.info paramsSummaryLog(workflow)
  
   // Write each accession string to its own file, named by sample id
   ch_ncbi = input_type.ncbi
      .collectFile { row -> [ "${row[0]}.txt", row[1] + '\n' ] }
      .map { f -> [ [id: f.baseName], f ] }

   NCBIGENOMEDOWNLOAD (
      ch_ncbi.map { meta, f -> meta },   // val meta
      ch_ncbi.map { meta, f -> f },      // path accessions (now an actual file)
      [],                                 // path taxids
      params.groups
   )

   // Build input channels
   ch_gff = NCBIGENOMEDOWNLOAD.out.gff
    .mix( input_type.local.map { [ [id: it[0]], file(it[2]) ] } )

   ch_fna_raw = NCBIGENOMEDOWNLOAD.out.fna
      .mix( input_type.local.map { [ [id: it[0]], file(it[1]) ] } )

   // Split on .gz, decompress only what needs it
   ch_fna_gz    = ch_fna_raw.filter { meta, fna -> fna.name.endsWith('.gz') }
   ch_fna_plain = ch_fna_raw.filter { meta, fna -> !fna.name.endsWith('.gz') }

   GUNZIP ( ch_fna_gz )

   ch_fna = GUNZIP.out.gunzip.mix( ch_fna_plain )

   // Convert GFF to AGAT format

   AGAT_CONVERTSPGXF2GXF ( ch_gff )

   // AGAT: gff channel + empty config
   AGAT_SPKEEPLONGESTISOFORM ( AGAT_CONVERTSPGXF2GXF.out.output_gff , [] )

   // Join fna + agat gff by meta, then split for GFFREAD's two inputs
   ch_fna_gff = ch_fna.join( AGAT_SPKEEPLONGESTISOFORM.out.gff )

   GFFREAD (
      ch_fna_gff.map { meta, fna, gff -> [ meta, gff ] },  // tuple val(meta), path(gff)
      ch_fna_gff.map { meta, fna, gff -> fna }             // path fasta (no meta)
   )

   FILTER_FASTA ( GFFREAD.out.gffread_fasta )

   // Use filtered fasta as merge_ch for all downstream processes
   merge_ch = FILTER_FASTA.out

   // If you wish to run GO, we use eggnogmapper:
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

      // Strip meta and collect plain go files for downstream use
      ch_go_files = EGGNOG_TO_GO.out.go_file
         .map { meta, go -> go }
         .collect()
   }

   if (params.stats){
      BUSCO_BUSCO (  GFFREAD.out.proteins_busco , 
                     params.busco_mode,
                     params.busco_lineage,
                     params.busco_lineages_path ?: [],
                     params.busco_config ?: [], 
                  )
      ch_versions = ch_versions.mix(BUSCO_BUSCO.out.versions.first())
   
      AGAT_SPSTATISTICS (  GFFREAD.out.gffs_agat  )
      ch_versions = ch_versions.mix(AGAT_SPSTATISTICS.out.versions.first())

      QUAST (  GFFREAD.out.fasta_quast,
               GFFREAD.out.gffs_agat
            )
      ch_versions = ch_versions.mix(QUAST.out.versions.first())
   }

   if (params.predownloaded_fasta && params.predownloaded_gofiles) {
      channel.fromPath(params.predownloaded_fasta)
         .mix( merge_ch.map { meta, fasta -> fasta } )
         .collect()
         .set{ proteins_ch }

      //ORTHOFINDER_GO ( proteins_ch.map { files -> [ [id: "ortho_go"], files ] }, [[],[]] )
   }
   else if ( params.ensembl_dataset && params.ensembl_biomart ){

      input_biomart = channel
         .value(params.ensembl_biomart)
         .ifEmpty { error "Cannot find the host name: ${params.ensembl_biomart}" }

      input_dataset = channel
         .fromPath(params.ensembl_dataset)
         .splitText().map{it -> it.trim()}
         .ifEmpty { error "Cannot find the dataset file: ${params.ensembl_dataset}" }

      GET_DATA ( input_biomart, input_dataset )
      ch_versions = ch_versions.mix(GET_DATA.out.versions)

      GET_DATA.out.gene_ontology_files.collect().set{ go_file_ch }

      GET_DATA.out.fasta_files
         .mix( merge_ch.map { meta, fasta -> fasta } )
         .collect()
         .set{ proteins_ch }

      //ORTHOFINDER_GO ( proteins_ch.map { files -> [ [id: "ortho_go"], files ] }, [[],[]] )
   }

   if (params.go_expansion) {
      //GO_EXPANSION ( GO_ASSIGN.out.go_counts.collect() )
      //ch_versions = ch_versions.mix(GO_EXPANSION.out.versions)
   }


   if (params.skip_cafe == null) {
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

      if (params.run_eggnog) {
        EGGNOG_TO_OG_GO (
            ch_go_files,                              // per-species go.txt files
            ORTHOFINDER_CAFE.out.orthologues          // Orthogroups.tsv
        )

        CAFE_GO (
            CAFE.out.result,
            CAFE.out.N0_table,
            EGGNOG_TO_OG_GO.out.og_go
        )
        ch_versions = ch_versions.mix(CAFE_GO.out.versions.first())
      }

   }


   if (params.chromo_go && params.run_eggnog) {
    CHROMO_GO (
        AGAT_SPKEEPLONGESTISOFORM.out.gff.map { meta, gff -> gff }.collect(),
        ch_go_files,
        ORTHOFINDER_CAFE.out.orthologues
    )
    ch_versions = ch_versions.mix(CHROMO_GO.out.versions)
   }


   CUSTOM_DUMPSOFTWAREVERSIONS ( ch_versions.collectFile(name: 'collated_versions.yml') )


}

workflow.onComplete { 
        println ( workflow.success ? "\nDone!\n" : "Oops... something went wrong" )
}

