#!/usr/bin/env nextflow
// Using DSL-2
nextflow.enable.dsl=2

process GENES_FROM_DISGENET {
  container "${params.container__r}"
  label 'low_proc'

  publishDir "${params.output_folder}/diseasemap", mode: 'copy', overwrite: true
                  
  input:
    path(disgenet)
    val(disease)
    val(threshold)

    
  output:
    path("geneList_disgenet.tsv")

  script:
    """
      get_disgenet_genes.R \
        ${disgenet} \
        ${disease} \
        ${threshold}

      echo "genes" > geneList_disgenet.tsv
      cat geneList.tsv >> geneList_disgenet.tsv
    """
}
