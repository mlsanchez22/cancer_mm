#!/usr/bin/env nextflow
// Using DSL-2
nextflow.enable.dsl=2

process CALCULATE_CIRCUIT_ACTIVITY {
  container "${params.container__r}"
  label 'low_proc'

  input:
    path(normalize_gene_counts)
    path(hipathia_rds)
    val(circuit_name) 

  output:    
    path("path_vals_transposed-*.rds")

  script:   
    """
      gene_expr_to_circ_activ.R \
        ${normalize_gene_counts} \
        ${hipathia_rds} \
        ${circuit_name}
    """
}
