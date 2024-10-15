#!/usr/bin/env nextflow

// Using DSL-2
nextflow.enable.dsl=2


/*
================================================================
	MODULES
================================================================
*/
//TODO: include modules
include {PARSE_DRUGBANK_XML} from "./modules/parse_drugbank_xml"
include {GENES_FROM_DISGENET} from "./modules/genes_from_disgenet"
include {PARSE_HIPATHIA} from "./modules/parse_hipathia"
include {ANNOTATE_HIPATHIA} from "./modules/annotate_hipathia"
include {
	GENE_ID_CONVERSION as HIPATHIA_CONV;
	GENE_ID_CONVERSION as DISGENET_CONV;
	GENE_ID_CONVERSION as DRUGBANK_CONV;
	GENE_ID_CONVERSION as GTEX_ENTREZ_CONV} from "./modules/gene_id_conversion"
include {GENES_TO_CIRCUITS} from "./modules/gene_to_circuits"
include {ANNOTATE_DRUGBANK} from "./modules/annotate_drugbank"
include {GET_GTEX_MAP} from "./modules/get_gtex_map"
include {NORMALIZE_EXPR as NORMALIZE_GTEX} from "./modules/normalize_expr"
include {CALCULATE_CIRCUIT_ACTIVITY as CIRCUIT_ACTIVITY_GTEX} from "./modules/calculate_circuit_activity"
include {CONCATENATE_CIRCUIT_ACTIVITIES} from "./modules/concatenate_circuit_activities"
include {CREATE_DREXML_ENV} from "./modules/create_exp_design"
include {RUN_DREXML} from "./modules/run_drexml"
include {
	COLLAPSE_COUNTS as COLLAPSE_GTEX_ENSEMBL;
	COLLAPSE_COUNTS as COLLAPSE_GTEX_SYMBOL} from "./modules/collapse_counts"
include {DRUGBANK_GENE_TARGETS} from "./modules/drugbank_gene_targets"
include {GTEX_SYMBOL_CONV} from "./modules/gtex_symbol_conv"
include {PLOT_DREXML} from "./modules/plot_drexml"


/*
================================================================
	SUBWORKFLOWS
================================================================
*/
//TODO: include subworkflows


/*
================================================================
	FUNCTIONS
================================================================
*/
// Function which prints help message text
//TODO: add help message
def helpMessage() {
    log.info"""
Usage:

nextflow run babelomics/schizorw <ARGUMENTS>
Either the disease map (gene list or circuit) or a disease map source (disgenet, differential activation, HPO, ...) is requiered.

	
	Drugbank preprocessing:
		--drugbank 				Drugbank compressed XML file as provided by https://go.drugbank.com/ (*.xml.gz)
		--atc					The Anatomical Therapeutic Chemical (ATC) Classification table from  https://bioportal.bioontology.org/ (*.csv)
		--drexml_target_genes	Target gene list file (one gene -in gene symbol format - per line, no header). If provided, drexml uses this list of target genes and ignores the drugbank and atc input parameters (skip drugbank preprocessing). Otherwise, the target gene list are retrieved from drugbank and ATC.
	
	GTEx preprocessing:
		--gtex					The Genotype-Tissue Expression (GTEx) RNA-Seq gene read counts from https://gtexportal.org/ (*.gct.gz)
		--drexml_gene_expr_norm		GTEx expression data normalized by TMM. Skip GTEx preprocessing (requires drexml_pathvals).
		--drexml_pathvals			GTEx circuits pathvals calculated by hipathia. Skip GTEx preprocessing (requires drexml_gene_expr_norm).

	Disease map circuits:
		--disease_map_circuits	Hipathia disease map circuit file (one circuit per line, no header). If provided, drexml uses this as the disease map source (overwrites others disease map sources).
		--disease_map_gene_list	Disease map gene list (one gene -in gene symbol format - per line, no header). If provided, the pipeline uses this file to retrieve the hipathia disease map circuits (overwrites others disease map sources).

	DISGENET disease map:
		--use_disgenet			Use DISGENET to retrieve the list of disease map circuits. Default TRUE. 
		--disgenet				Disgenet database as provided by https://www.disgenet.org/ (*.tsv.gz)
		--disgenet_disease_id	Disgenet disease id. This parameter is used to filter out disease map genes from disgenet.
		--disgenet_threshold	Percentile in range 0-1 above which select highest scoring genes. Default 0.75

	Differential activation disease map:
		--use_diff_act			Performs a case/control circuit differential activation test to retrieve the list of disease map circuits. ***UNDER DEVELOPMENT***

	HPO disease map:
		--use_hpo			Use the HPO database to retrieve the list of disease map circuits. ***UNDER DEVELOPMENT***
	
	Computation params:		
		--n_gpus <int>				Use the indicated number of GPUS (0 for CPUs). ***UNDER DEVELOPMENT***

	Output:
		--output_folder				Results output folder (Default: ${launchDir}/results/)
""".stripIndent()
}


/*
================================================================
DEFINE SUBWORKFLOW DEPENDENCES
================================================================
*/
DRUGBANK_PREPROC_DEPENDENCES = [
	params.drugbank,
	params.atc
]
GTEX_PREPROC_DEPENDENCES = [
	params.gtex
]
DISGENET_PREPROC_DEPENDENCES = [
	params.disgenet,
	params.disgenet_disease_id,
	params.disgenet_threshold
]

/*
================================================================
MAIN WORKFLOW
================================================================
*/
workflow {
	//**********
	//*  INIT  *
	//**********
	// Print input parameters
	log.info "Input parameters:"
	params.each { key, value ->
		if (value) {println "${key}: ${value}"} }
	log.info "--------------------------\n"


	// set channels
	ch_disease_map_circuits		= params.disease_map_circuits 	? channel.fromPath(params.disease_map_circuits)		: channel.empty()
	ch_disease_map_gene_list	= params.disease_map_gene_list	? channel.fromPath(params.disease_map_gene_list)	: channel.empty()
	ch_drexml_target_genes		= params.drexml_target_genes 	? channel.fromPath(params.drexml_target_genes) 		: channel.empty()
	ch_drexml_gene_expr_norm	= params.drexml_gene_expr_norm 	? channel.fromPath(params.drexml_gene_expr_norm) 	: channel.empty()
	ch_drexml_pathvals			= params.drexml_pathvals 		? channel.fromPath(params.drexml_pathvals) 			: channel.empty()
	ch_drexml_results			= params.drexml_results 		? channel.fromPath(params.drexml_results) 			: channel.empty()
	ch_pathways					= channel.empty()

	// define disease map building method
	use_disgenet	= (params.disease_map_circuits || params.disease_map_gene_list ? false : params.use_disgenet)
	use_diff_act 	= false //FIXME: under development
	use_hpo			= false //FIXME: under development

	// set preprocessing steps
	preproc_drugbank	= !params.drexml_target_genes ? true : false
	preproc_hipathia	= !params.drexml_gene_expr_norm 
							|| !params.drexml_pathvals 
							|| !params.disease_map_circuits 
							|| !params.disease_map_gene_list ? true : false
	preproc_gtex		= !params.drexml_gene_expr_norm 
							|| !params.drexml_pathvals ? true : false
	preproc_drexml		= !params.drexml_results ? true : false

	
	// Show help message if the user specifies the --help flag at runtime
	if (params.help){
		helpMessage()
		exit 1
	}

	//**************************
	//* CHECK INPUT PARAMETERS *
	//**************************
	// configure the disease map source
	if (use_disgenet && DISGENET_PREPROC_DEPENDENCES.any { it == false }){
		log.error "ERROR: some DISGENET preprocessing dependences coun't be found!."
		helpMessage()
		exit 1
	}	

	// Check drugbank preprocessing dependences
	if (preproc_drugbank && DRUGBANK_PREPROC_DEPENDENCES.any { it == false }){
		log.error "ERROR: Either the gene target list or the drugbank preprocessing dependences are needed."
		helpMessage()
		exit 1
	}

	// Check GTEx preprocessing dependences
	if (preproc_gtex && GTEX_PREPROC_DEPENDENCES.any { it == false }){
		log.error "ERROR: Either GTEx normalized expression matrix and circuit activity values or GTEx preprocessing dependences are needed."
		helpMessage()
		exit 1
	}

	// Check disease map source
	if (!params.disease_map_circuits && !params.disease_map_gene_list && !use_disgenet){
		log.error "ERROR: Either the disease map (list of circuits or genes) or a disease map building method is needed."
		helpMessage()
		exit 1
	}



	//**************************
	//* DRUGBANK PREPROCESSING *
	//**************************
	if(preproc_drugbank){
		// parse drugbank xml file
		PARSE_DRUGBANK_XML(
			params.drugbank
		) 

		// translate uniprot_id to gene symbol
		DRUGBANK_CONV(
			PARSE_DRUGBANK_XML.out,
			"uniprot_id",
			"gene_symbol",
			"uniprot_id",
			"gene_symbol"
		)

		// add ATC codes and simplified drugs actions
		ANNOTATE_DRUGBANK(
			DRUGBANK_CONV.out,
			params.atc
		)
			
		ch_drexml_target_genes = DRUGBANK_GENE_TARGETS(
			ANNOTATE_DRUGBANK.out
		)

	} else {
		log.info "INFO: target genes taken from ${params.drexml_target_genes}"
	}


	//**************************
	//* HIPATHIA PREPROCESSING *
	//**************************
	// Preprocess hipathia if needed
	if(preproc_hipathia){
		// Extract hipathia metaginfo
		PARSE_HIPATHIA()

		// Translate entrez_id to gene symbol
		HIPATHIA_CONV(		
			PARSE_HIPATHIA.out.hipathia_tsv,
			"gene",
			"gene_symbol",
			"entrez_id",
			"gene_symbol"
		) 

		// add some other usefull columns
		ANNOTATE_HIPATHIA(
			HIPATHIA_CONV.out
		)

		hipathia_annot 			= HIPATHIA_CONV.out
		hipathia_metaginfo_rds 	= PARSE_HIPATHIA.out.hipathia_rds
	}	




	//**********************************
	//* DISEASE CIRCUITS PREPROCESSING *
	//**********************************	
	// Retrieve the list of disease associated hipathia circuits	
	if(!params.disease_map_circuits){
		if(!params.disease_map_gene_list){
			if(use_disgenet){
				log.info "INFO: Using DISGENET to build the disease map"

				//******** DISGENET DISEASE ASSOCIATED GENE LIST ********
				// Convert entrez_id to gene_symbol (just in case)
				DISGENET_CONV(			
					params.disgenet,
					"geneId",
					"gene_symbol",
					"entrez_id",
					"gene_symbol"
				) 
				
				// retrieve the list of disease associated genes
				GENES_FROM_DISGENET(
					DISGENET_CONV.out,
					params.disgenet_disease_id,
					params.disgenet_threshold
				)
						
				ch_disease_map_gene_list = GENES_FROM_DISGENET.out

			} else if(use_diff_act){
				log.error "ERROR: Differential activation disease map construction is not implemented yet!"
				exit 1

			} else if(use_hpo){
				log.error "ERROR: HPO disease map construction is not implemented yet!"
				exit 1

			} else {
				log.error "ERROR: At least a disease map building method must be specified!"
				exit 1
			}

		} else {								
			log.info "INFO: Disease map will be built from the gene list ${params.disease_map_gene_list}"	
		}

		// Retrieve the list of hipathia circuits from a set of genes. 
		GENES_TO_CIRCUITS(		
			hipathia_annot,
			ch_disease_map_gene_list,
			"gene_symbol",	
			"circuit"
		)

		ch_disease_map_circuits = GENES_TO_CIRCUITS.out.circuits_drexml
			

	} else {		
		log.info "INFO: Disease map circuit list taken from ${params.disease_map_circuits}"
	}




	//**********************
	//* GTEx PREPROCESSING *
	//**********************
	if (preproc_gtex){
		// Collapse ESNG from GTEx (just in case) and remove the first 2 lines
		COLLAPSE_GTEX_ENSEMBL(
			params.gtex,
			"Name",
			2)
		
		// Translate ensembl ids to gene symbols
		GTEX_ENTREZ_CONV(
			COLLAPSE_GTEX_ENSEMBL.out,
			"Name",
			"entrez_id",
			"ensembl_id",
			"entrez_id"
		)
		
		// Collapse genes counts (add counts)
		COLLAPSE_GTEX_SYMBOL(
			GTEX_ENTREZ_CONV.out,
			"entrez_id",
			0)
	
		// Normalize gtex gene counts (TMM + cpm)
		//TODO: add rds to speed up circuit activity calculation
		NORMALIZE_GTEX(
			COLLAPSE_GTEX_SYMBOL.out
		)

		// create channel to compute disease pathways activities
		ch_pathways = ch_disease_map_circuits
				.splitCsv(header: true, strip: true, sep: '\t')			
				.map { row -> row['index'].find("hsa[0-9]+") }						
				.unique()
				.dump(tag: 'processed pathways: ')


		// Calculate pathways activity in parallel
		CIRCUIT_ACTIVITY_GTEX(
			NORMALIZE_GTEX.out.gene_expr_x_samples_rds,
			hipathia_metaginfo_rds,
			ch_pathways
		)

		// Combine results
		CONCATENATE_CIRCUIT_ACTIVITIES(
			CIRCUIT_ACTIVITY_GTEX.out.collect()
		)

		// Translate (normalized gene expression) entrez_ids to gene_symbols
		GTEX_SYMBOL_CONV(
			NORMALIZE_GTEX.out.samples_x_gene_expr
		)

		ch_drexml_gene_expr_norm	= GTEX_SYMBOL_CONV.out
		ch_drexml_pathvals			= CONCATENATE_CIRCUIT_ACTIVITIES.out.samples_x_pathvals

	} else {
		log.info "INFO: GTEx expression matrix taken from ${params.drexml_gene_expr_norm}"
		log.info "INFO: GTEx circuit activity taken from ${params.drexml_pathvals}"
	}




	//********************
	//* DREXML EXECUTION *
	//********************
	if (!params.drexml_results){
		// create drexml environment file
		CREATE_DREXML_ENV(
			ch_drexml_gene_expr_norm,	
			ch_drexml_pathvals,
			ch_disease_map_circuits,
			ch_drexml_target_genes		 
		)

		// run drexml
		RUN_DREXML(
			ch_drexml_gene_expr_norm,
			ch_drexml_pathvals,
			ch_disease_map_circuits,
			ch_drexml_target_genes,
			CREATE_DREXML_ENV.out
		)	

		ch_drexml_results = RUN_DREXML.out

	} else {
		log.info "INFO: Drexml results taken from ${params.drexml_results}"
	}

	PLOT_DREXML(
		ch_drexml_results
	)
}