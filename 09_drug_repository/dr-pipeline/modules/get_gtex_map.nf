#!/usr/bin/env nextflow
// Using DSL-2
nextflow.enable.dsl=2

process GET_GTEX_MAP {
  container "${params.container__r}"
  label 'low_proc'

  //publishDir "${params.output_folder}", mode: 'copy', overwrite: true

  input:
    path(gtex)
    var(column_from)
    var(column_to)
    var(format_from)
    var(format_to)

  output:
    path("gtex_map.tsv") 

  script:   
    """
      gene_id_conversion.R \
        ${gtex} \
        ${column_from} \
        ${column_to} \
        ${format_from} \
        ${format_to} |  awk \
          NR==1 {
            for (i=1; i<=NF; i++) {
              ix[$i] = i
            }
            print "${column_from}","${column_to}"
          }
          
          NR>1 {
            print $ix["${column_from}"],$ix["${column_to}"]
          } OFS="\t" > gtex_map.tsv
    """
}
