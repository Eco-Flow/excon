#!/usr/bin/env nextflow

log.info """\
 =========================================

 EXCON (v3.0)

 -----------------------------------------

 Authors:
   - Chris Wyatt <c.wyatt@ucl.ac.uk>
   - Simon Murray <simon.murray@ucl.ac.uk>

 -----------------------------------------

 Copyright (c) 2021

 =========================================""".stripIndent()

include { GET_DATA } from './modules/local/getdata.nf'
include { ORTHOFINDER as ORTHOFINDER_GO} from './modules/local/orthofinder.nf'
include { ORTHOFINDER as ORTHOFINDER_CAFE } from './modules/local/orthofinder.nf'
include { GO_ASSIGN } from './modules/local/go_assign.nf'
include { GO_EXPANSION  } from './modules/local/go_expansion.nf'
include { DOWNLOAD_NCBI } from './modules/local/download_ncbi.nf'
include { GFFREAD } from './modules/local/gffread.nf'
include { CAFE } from './modules/local/cafe.nf'
include { CHROMO_GO } from './modules/local/chromo_go.nf'
include { CAFE_GO } from './modules/local/cafe_go.nf'

include { validateParameters; paramsHelp; paramsSummaryLog } from 'plugin/nf-validation'
include { CUSTOM_DUMPSOFTWAREVERSIONS } from './modules/nf-core/custom/dumpsoftwareversions/main'


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

   DOWNLOAD_NCBI ( input_type.ncbi )

   ch_versions = ch_versions.mix(DOWNLOAD_NCBI.out.versions.first())

   GFFREAD ( DOWNLOAD_NCBI.out.genome.mix(input_type.local) )

   ch_versions = ch_versions.mix(GFFREAD.out.versions.first())

   merge_ch = GFFREAD.out.longest.collect()
   
   // If you have precomiled GO data, or want to download from biomart else do not do GO enrichmetn analysis:
   if (params.predownloaded_fasta && params.predownloaded_gofiles) {
      channel.fromPath(params.predownloaded_fasta).mix(merge_ch).collect().set{ proteins_ch }
      channel.fromPath(params.predownloaded_gofiles).collect().set{ go_file_ch }
      ORTHOFINDER_GO ( proteins_ch )
      ch_versions = ch_versions.mix(ORTHOFINDER_GO.out.versions)

      if (params.predownloaded_fasta || params.ensembl_dataset){
         GO_ASSIGN ( go_file_ch , ORTHOFINDER_GO.out.orthologues, GFFREAD.out.longest , GFFREAD.out.gene_to_isoforms.collect() )
         //ch_versions = ch_versions.mix(GO_ASSIGN.out.versions.first())
      }
   }
   else if (params.ensembl_dataset && params.ensembl_host && params.ensembl_repo){
      background_species = channel
         .fromPath(params.ensembl_dataset) 
         .splitText().map{it -> it.trim()}
         .ifEmpty { error "Cannot find the dataset file: ${params.ensembl_dataset}" }

      input_host = channel
         .value(params.ensembl_host)
         .ifEmpty { error "Cannot find the host name: ${params.ensembl_host}" }

      input_repo = channel
         .value(params.ensembl_repo)
         .ifEmpty { error "Cannot find the repo name: ${params.ensembl_repo}" }

      GET_DATA ( input_repo, input_host, background_species )
      ch_versions = ch_versions.mix(GET_DATA.out.versions)

      GET_DATA.out.gene_ontology_files.set{ go_file_ch }

      GET_DATA.out.fasta_files.mix(merge_ch).collect().set{ proteins_ch }

      ORTHOFINDER_GO ( proteins_ch )
      ch_versions = ch_versions.mix(ORTHOFINDER_GO.out.versions)

      if (params.predownloaded_fasta || params.ensembl_dataset){
         GO_ASSIGN ( go_file_ch , ORTHOFINDER_GO.out.orthologues, GFFREAD.out.longest , GFFREAD.out.gene_to_isoforms.collect() )
         //ch_versions = ch_versions.mix(GO_ASSIGN.out.versions.first())
      }
   }
   else {
      //User chooses not to run GO enrichment analysis
   }

   //Run GO expansion analysis
   if (params.go_expansion) {
     GO_EXPANSION ( GO_ASSIGN.out.go_counts.collect() )
     ch_versions = ch_versions.mix(GO_EXPANSION.out.versions)
   }

   //Run chromosome GO analysis
   if (params.chromo_go) {
     CHROMO_GO ( GFFREAD.out.gffs.collect() , GO_ASSIGN.out.go_hash.collect() , ORTHOFINDER.out.orthologues )
     ch_versions = ch_versions.mix(CHROMO_GO.out.versions)
   }

   //Run Orthofinder for CAFE using just input (focal) species.
   ORTHOFINDER_CAFE ( merge_ch )
   //No need to collect versions from orthofinder module twice

   //Run Cafe analysis of expanded and contracted gene families.
   CAFE ( ORTHOFINDER_CAFE.out.no_ortho, ORTHOFINDER_CAFE.out.speciestree )
   ch_versions = ch_versions.mix(CAFE.out.versions)

   if (params.predownloaded_fasta || params.ensembl_dataset) {
      CAFE_GO ( CAFE.out.result, CAFE.out.N0_table, GO_ASSIGN.out.go_og )
      ch_versions = ch_versions.mix(CAFE_GO.out.versions.first())
   }
   else{
      //Do nothing no go data to run this part. 
   }

   CUSTOM_DUMPSOFTWAREVERSIONS ( ch_versions.collectFile(name: 'collated_versions.yml') )

}