#!/usr/bin/env nextflow

// Using DSL-2
nextflow.enable.dsl=2

process PARSE_DRUGBANK_XML {
  container "${params.container__python}"
  label 'low_proc'

  //publishDir "${params.output_folder}", mode: 'copy', overwrite: true
                  
  input:
    path(drugbank_xml_file)

  output:
    path "drugbank.tsv"

  script:   
    """
      parser.py parse ${drugbank_xml_file} drugbank.tsv
    """
}
