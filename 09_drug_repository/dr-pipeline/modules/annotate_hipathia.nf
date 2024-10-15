#!/usr/bin/env nextflow
// Using DSL-2
nextflow.enable.dsl=2

process ANNOTATE_HIPATHIA {
  container "${params.container__r}"
  label 'low_proc'

  publishDir "${params.output_folder}/preproc/hipathia", mode: 'copy', overwrite: true, saveAs: {"hipathia_circuits_annotated.tsv"} 

  input:
    path(hipathia_circuits) 

  output:
    path("*-annot.tsv")

  script:   
    """
      get_hipathia_genes_annot2.R \
        ${hipathia_circuits}
    """
}
