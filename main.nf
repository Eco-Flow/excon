/*
 * Copyright (c) 2021
 */


/*
 * Authors:
 * - Chris Wyatt <chris.wyatt@seqera.io>
 */

/*
 * Default pipeline parameters (on test data). They can be overriden on the command line eg.
 * given `params.input` specify on the run command line `--input /path/to/input`.
 */

params.input= "data/Example.csv"
params.ensembl_repo="metazoa_mart"
params.ensembl_host='https://metazoa.ensembl.org'
params.ensembl_dataset="example.txt"
params.predownloaded_fasta= "./Background_species_folder/*"
params.predownloaded_gofiles= "./Background_gofiles_folder/*"
params.outdir = "results"
params.download= false

params.cafe= false
params.chromo_go= false
params.orthofinder= false
params.go_assign= false
params.go_expansion= false
params.cafe_go= false

//For CPU and Memory of each process: see conf/docker.config

log.info """\
 ===================================
	GOATEE v2.0
 ===================================
 input file                           : ${params.input}
 list of background species           : ${params.ensembl_dataset}
 out directory                        : ${params.outdir}
 """

//================================================================================
// Include modules
//================================================================================

include { GET_DATA } from './modules/getdata.nf'
include { ORTHOFINDER } from './modules/orthofinder.nf'
include { ORTHOFINDER as ORTHOFINDER_2 } from './modules/orthofinder.nf'
include { GO_ASSIGN } from './modules/go_assign.nf'
include { GO_EXPANSION  } from './modules/go_expansion.nf'
include { DOWNLOAD_NCBI } from './modules/download_ncbi.nf'
include { GFFREAD } from './modules/gffread.nf'
include { CAFE } from './modules/cafe.nf'
include { CHROMO_GO } from './modules/chromo_go.nf'
include { CAFE_GO } from './modules/cafe_go.nf'



Channel
    .fromPath(params.input)
    .splitCsv()
    .branch { 
        ncbi: it.size() == 2 
        local: it.size() == 3
    }
    .set { input_type }

workflow {

	DOWNLOAD_NCBI ( input_type.ncbi )

	GFFREAD ( DOWNLOAD_NCBI.out.genome.mix(input_type.local) )

	merge_ch = GFFREAD.out.longest.collect()

	if (params.go_assign || params.cafe_go ){
		if (params.download){
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

			GET_DATA.out.gene_ontology_files.set{ go_file_ch }

			GET_DATA.out.fasta_files.mix(merge_ch).collect().set{ proteins_ch }

		}
		else{
			channel.fromPath(params.predownloaded_fasta).mix(merge_ch).collect().set{ proteins_ch }
			channel.fromPath(params.predownloaded_gofiles).collect().set{ go_file_ch }
		}

		ORTHOFINDER ( proteins_ch )

		GO_ASSIGN ( go_file_ch , ORTHOFINDER.out.orthologues, GFFREAD.out.longest , GFFREAD.out.gene_to_isoforms.collect() )

		if (params.go_expansion){
			GO_EXPANSION ( GO_ASSIGN.out.go_counts.collect() )
		}

		if (params.chromo_go){
			CHROMO_GO ( GFFREAD.out.gffs.collect() , GO_ASSIGN.out.go_hash.collect() , ORTHOFINDER.out.orthologues )
		}


	}

	if (params.orthofinder){

		ORTHOFINDER_2 ( merge_ch )

		if (params.cafe){

			CAFE ( ORTHOFINDER_2.out.no_ortho  , ORTHOFINDER_2.out.speciestree )	

			if (params.cafe_go){
				CAFE_GO ( CAFE.out.result , CAFE.out.N0_table , GO_ASSIGN.out.go_og.first() )
			}	
		}

	}
	else{

		if (params.cafe){

			ORTHOFINDER_2 ( merge_ch )
			CAFE ( ORTHOFINDER_2.out.no_ortho  , ORTHOFINDER_2.out.speciestree )	

			if (params.cafe_go){
                                CAFE_GO ( CAFE.out.result , CAFE.out.N0_table , GO_ASSIGN.out.go_og.first() )
                        }
		}
	}

}

workflow.onComplete {
	println ( workflow.success ? "\nDone! Open your report in your browser --> $params.outdir/report.html (if you added -with-report flag)\n" : "Hmmm .. something went wrong" )
}

