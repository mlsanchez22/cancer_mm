#!/usr/bin/env nextflow
// Using DSL-2
nextflow.enable.dsl=2

process GTEX_SYMBOL_CONV {
  container "${params.container__python}"
  label 'low_proc'

  publishDir "${params.output_folder}/preproc/GTEx", mode: 'copy', overwrite: true, saveAs: {"GTEx_norm_symbol.tsv"}

  input:
    path(samples_x_genes)

  output:
    path("samples_x_gene_symbol.tsv")

  script:   
    """
      # Transform the header in one-colum file
      for i in \$(head -n 1 ${samples_x_genes})  
      do
        echo \$i
      done > gene_list_entrez.tsv

      # convert entrez_ids to gene_symbol
      gene_id_conversor.py translate \
        --skip_missing true \
        --act_on_duplicated duplicate_entry \
        --in_col_name index \
        --out_col_name gene_symbol \
        --in_id_format entrez_id \
        --out_id_format gene_symbol \
        gene_list_entrez.tsv \
        gene_list_symbol.tsv

      # Create the new header
      a="index"
      for i in \$(tail -n+2 gene_list_symbol.tsv|cut -f 2)
      do
        a="\$a\t\$i"
      done

      echo -e "\${a}" > samples_x_gene_symbol.tsv

      # append the old matrix
      tail -n+2 ${samples_x_genes} >> samples_x_gene_symbol.tsv
    """
}
