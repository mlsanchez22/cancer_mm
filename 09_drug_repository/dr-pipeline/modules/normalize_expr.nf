#!/usr/bin/env nextflow
// Using DSL-2
nextflow.enable.dsl=2

process NORMALIZE_EXPR {
  container "${params.container__r}"
  label 'low_proc'

  //publishDir "${params.output_folder}", mode: 'copy', overwrite: true

  input:
    path(gene_counts) 

  output:
    path("tmm-normalized_expr.rds", emit: "gene_expr_x_samples_rds")
    path("tmm-normalized_expr_transposed.tsv", emit: "samples_x_gene_expr")

  script:   
    """
      # gene_expr_norm <counts> <low-expression gene filtering> <log2 transformation of gene expr>
      gene_expr_norm.R \
        ${gene_counts} \
        "TRUE" \
        "TRUE"
    """
}
