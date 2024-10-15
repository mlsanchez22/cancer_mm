#!/usr/bin/env nextflow
// Using DSL-2
nextflow.enable.dsl=2

process GENE_ID_CONVERSION {
  container "${params.container__python}"
  label 'high_mem'

  //publishDir "${params.output_folder}", mode: 'copy', overwrite: true

  input:
    path(input_file) 
    val(column_from)
    val(column_to)
    val(format_from)
    val(format_to)    

  output:
    path("annotated_gene_matrix.tsv")

  script:   
    """
      gene_id_conversor.py translate \
        --skip_missing true \
        --act_on_duplicated duplicate_entry \
        --in_col_name ${column_from} \
        --out_col_name ${column_to} \
        --in_id_format ${format_from} \
        --out_id_format ${format_to} \
        ${input_file} \
        annotated_gene_matrix.tsv
    """
}
