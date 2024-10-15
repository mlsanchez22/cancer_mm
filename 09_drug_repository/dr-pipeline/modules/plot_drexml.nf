#!/usr/bin/env nextflow

// Using DSL-2
nextflow.enable.dsl=2

process PLOT_DREXML {
  container "${params.container__drexml}"
  label 'low_proc'  //FIXME

  publishDir "${params.output_folder}/drexml/plots", mode: 'copy', overwrite: true
                  
  input:
    path(drexml_results)

  output:
    path("*.pdf")

  script:   
    """
      export PYSTOW_NAME="`pwd`/.data"
      export MPLCONFIGDIR="`pwd`/.data"
      
      drexml plot \
        ${drexml_results}/shap_selection_symbol.tsv \
        ${drexml_results}/shap_summary_symbol.tsv \
        ${drexml_results}/stability_results_symbol.tsv \
        .
    """
}
