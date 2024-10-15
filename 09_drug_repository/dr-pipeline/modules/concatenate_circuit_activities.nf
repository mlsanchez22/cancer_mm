#!/usr/bin/env nextflow
// Using DSL-2
nextflow.enable.dsl=2

process CONCATENATE_CIRCUIT_ACTIVITIES {
  container "${params.container__r}"
  label 'single_proc'

  //FIXME: take care of file overwritting when this module is used to retrieve the disease map circuits
  publishDir "${params.output_folder}/preproc/GTEx", mode: 'copy', overwrite: true

  input:
    path(circuit_activity_files)

  output:    
    path("path_vals_transposed_concatenated.tsv", emit: "samples_x_pathvals")

  //TODO: add script name
  script:   
    """
      conca_path_vals_tt_files.R ${circuit_activity_files.join(' ')}
    """
}
