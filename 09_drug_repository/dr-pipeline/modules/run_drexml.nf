#!/usr/bin/env nextflow

// Using DSL-2
nextflow.enable.dsl=2

process RUN_DREXML {
  container "${params.container__drexml}"
  label 'high_proc'  

  publishDir "${params.output_folder}/drexml", 
    mode: 'copy', 
    overwrite: true, 
    pattern: 'results', 
    saveAs: { "tables" }
                  
  input:
    path(hipathia_norm)
    path(hipathia_pathvals)
    path(disease_circuits)
    path(drugbank_targets)
    path(drexml_env)

  output:
    path("results", emit: "results")

  script:   
    """
      export PYSTOW_NAME="`pwd`/.data"
      export MPLCONFIGDIR="`pwd`/.data"

      # FIXME: use generated dictionaries and resources
      # FIXME: use GPU docker container if n_gpus > 0
      drexml run \
        --n-gpus ${params.n_gpus} \
        --n-cpus -1 \
        ${drexml_env}

      rm -rf results/tmp
    """
}
