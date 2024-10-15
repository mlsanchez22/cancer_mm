#!/usr/bin/env nextflow
// Using DSL-2
nextflow.enable.dsl=2

process PARSE_HIPATHIA {
  container "${params.container__r}"
  label 'low_proc'

  publishDir "${params.output_folder}/preproc/hipathia", mode: 'copy', overwrite: true, pattern: "*.rds"
    
  output:
    path("pathways-hipathia*.rds", emit: "hipathia_rds")
    path("circuitGenes_hipathia*.tsv", emit: "hipathia_tsv")

  script:   
    """
      export ANNOTATION_HUB_CACHE="ahub_cache"

      get_hipathia_genes_annot1.R
    """
}
