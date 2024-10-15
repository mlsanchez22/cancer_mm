#!/usr/bin/env nextflow
// Using DSL-2
nextflow.enable.dsl=2

process ANNOTATE_DRUGBANK {
  container "${params.container__r}"
  label 'low_proc'

  publishDir "${params.output_folder}/preproc/drugbank", mode: 'copy', overwrite: true, saveAs: {"drugbank_simpl_atc.tsv"} 

  input:
    path(drugbank)
    path(atc_annots)

  output:
    path("drugbank_db-simpl_atc.tsv")

  script:   
    """
      drugbank_annotation.R \
      ${drugbank} \
      ${atc_annots} 
    """
}