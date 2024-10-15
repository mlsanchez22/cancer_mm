#!/usr/bin/env nextflow

// Using DSL-2
nextflow.enable.dsl=2

process CREATE_DREXML_ENV {
  label 'low_proc'

  //publishDir "${params.output_folder}/drexml", mode: 'copy', overwrite: true
                  
  input:
    path(hipathia_norm)
    path(hipathia_pathvals)
    path(disease_circuits)
    path(drugbank_targets)

  output:
    path "exp_design.env"

  shell:   
    """
      cat <<EOF > exp_design.env
###### EXPERIMENT DESIGN  ######
data_path="./"
gene_exp=!{hipathia_norm}
pathvals=!{hipathia_pathvals}
circuits=!{disease_circuits}
circuits_column="in_disease"
genes=!{drugbank_targets}
genes_column="is_target"
EOF
    """
}
