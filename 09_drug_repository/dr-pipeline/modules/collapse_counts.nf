#!/usr/bin/env nextflow

// Using DSL-2
nextflow.enable.dsl=2

process COLLAPSE_COUNTS {
  container "${params.container__python}"
  label 'high_mem'

  //publishDir "${params.output_folder}", mode: 'copy', overwrite: true
                  
  input:
    path(input_file)
    val(column)
    val(skiprows)

  output:
    path "counts.tsv"

  script:   
    """
      collapse.py \
        --how sum \
        --id_column ${column} \
        --skiprows ${skiprows} \
        ${input_file} \
        counts.tsv
    """
}
