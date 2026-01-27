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

include { GET_DATA } from './modules/local/getdata.nf'
include { ORTHOFINDER as ORTHOFINDER_GO} from './modules/nf-core/orthofinder/main.nf'
include { ORTHOFINDER as ORTHOFINDER_CAFE } from './modules/nf-core/orthofinder/main.nf'
include { GO_ASSIGN } from './modules/local/go_assign.nf'
include { GO_EXPANSION  } from './modules/local/go_expansion.nf'
include { NCBIGENOMEDOWNLOAD } from './modules/nf-core/ncbigenomedownload/main.nf'
include { GFFREAD } from './modules/local/gffread.nf'
include { CAFE } from './modules/local/cafe_with_retry.nf'
include { RESCALE_TREE } from './modules/local/rescale_tree.nf'
include { CHROMO_GO } from './modules/local/chromo_go.nf'
include { CAFE_GO } from './modules/local/cafe_go.nf'
include { CAFE_PLOT } from './modules/local/cafe_plot.nf'

include { validateParameters; paramsHelp; paramsSummaryLog } from 'plugin/nf-validation'
include { CUSTOM_DUMPSOFTWAREVERSIONS } from './modules/nf-core/custom/dumpsoftwareversions/main.nf'
include { BUSCO_BUSCO } from './modules/nf-core/busco/busco/main.nf'
include { AGAT_SPSTATISTICS } from './modules/nf-core/agat/spstatistics/main.nf'
include { QUAST } from './modules/nf-core/quast/main.nf'

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
  
   NCBIGENOMEDOWNLOAD ( input_type.ncbi.map { it[0] }, input_type.ncbi.map { it[1] }, [], params.groups)
   ch_versions = ch_versions.mix(NCBIGENOMEDOWNLOAD.out.versions.first())  

   GFFREAD ( NCBIGENOMEDOWNLOAD.out.fna.mix( input_type.local.map { [it[0],file(it[1])] } ), NCBIGENOMEDOWNLOAD.out.gff.mix(input_type.local.map { [it[0],file(it[2])] 
} ) ) 
   ch_versions = ch_versions.mix(GFFREAD.out.versions.first())

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
               //reference
               GFFREAD.out.gffs_agat
            )
      ch_versions = ch_versions.mix(QUAST.out.versions.first())
   }

   merge_ch = GFFREAD.out.longest.collect()
   
   
   if (params.skip_cafe == null) {
      //Run Orthofinder for CAFE using just input (focal) species.
      ORTHOFINDER_CAFE ( merge_ch.map { ["ortho_cafe", it] } )
      //No need to collect versions from orthofinder module twice

      //Rescale tree branch lengths to prevent numerical precision issues
      RESCALE_TREE ( ORTHOFINDER_CAFE.out.speciestree )

      //Run Cafe analysis of expanded and contracted gene families.
      CAFE ( ORTHOFINDER_CAFE.out.no_ortho, RESCALE_TREE.out.rescaled_tree )
      ch_versions = ch_versions.mix(CAFE.out.versions)

      CAFE_PLOT ( CAFE.out.result )
      ch_versions = ch_versions.mix(CAFE_PLOT.out.versions)

      if (params.predownloaded_fasta || params.ensembl_dataset) {
         CAFE_GO ( CAFE.out.result, CAFE.out.N0_table, GO_ASSIGN.out.go_og )
         ch_versions = ch_versions.mix(CAFE_GO.out.versions.first())
      }
   }

   CUSTOM_DUMPSOFTWAREVERSIONS ( ch_versions.collectFile(name: 'collated_versions.yml') )
}

workflow.onComplete { 
        println ( workflow.success ? "\nDone!\n" : "Oops... something went wrong" )
}

