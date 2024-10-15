#!/usr/bin/env nextflow
// Using DSL-2
nextflow.enable.dsl=2

process GENES_TO_CIRCUITS {
  container "${params.container__r}"
  label 'low_proc'

  publishDir "${params.output_folder}/diseasemap", mode: 'copy', overwrite: true
  
  input:
    path(gene_circuit_table)
    path(gene_list)
    val(gene_column_name)
    val(circuit_column_name)

  output:
    path("circuits_drexml.tsv", emit: "circuits_drexml")
    path("pathways.tsv", emit: "pathways")

  script:   
    """
      get_circuits_from_genes.R \
        ${gene_circuit_table} \
        ${gene_list} \
        ${gene_column_name} \
        ${circuit_column_name} 
             
      awk -F "\t" 'BEGIN{print "index\tin_disease"} {print \$1"\tTRUE"}' circuitList.tsv > circuits_drexml.tsv

      egrep -o 'hsa[0-9]+' circuitList.tsv |sort|uniq > pathways.tsv
    """
}

