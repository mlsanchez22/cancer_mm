#!/usr/bin/env nextflow
// Using DSL-2
nextflow.enable.dsl=2

process DRUGBANK_GENE_TARGETS {
  container "${params.container__r}"
  label 'low_proc'

  publishDir "${params.output_folder}/preproc/drugbank", mode: 'copy', overwrite: true, saveAs: {"drugbank_gene_targets.tsv"} 

  input:
    path(drugbank_annot)

  output:
    path("drugbank_targets.tsv")

  script:   
    """
      awk -F "\t" 'NR>1{print \$17"\tTRUE"}' ${drugbank_annot} \
        |sort \
        |uniq \
        |awk -F "\t" 'BEGIN{print "index\tis_target"} {print}' > drugbank_targets.tsv
    """
}