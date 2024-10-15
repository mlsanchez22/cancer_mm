#!/usr/bin/env python3
# coding: utf-8
"""
gene_id conversion.
Adapted from https://github.com/babelomics/dr-pipeline/blob/master/modules/resources/usr/bin/parser.py
"""

"""
Input Parameters:

    tsv_file_path: The TSV file that needs to be transformed.
    in_col_name: The name of the column containing the identifiers to be transformed. (NB: duplicated genes are not allowed)
    out_col_name: The name of the column containing the converted identifiers.
    in_id_format: The format of the input identifier [ensembl_id/uniprot_id/entrez_id].
    out_id_format: The format of the output identifier [entrez_id/gene_symbol].

Output:

    tsv_conv_file_path: A TSV file similar to the input file but with a new column containing the converted identifier.
"""
from pathlib import Path
import click
import pandas as pd
from biothings_client import get_client
import warnings
import sys

THIS_VERSION = 0.1 #not yet finished  :S

def read_file_df(tsv_file_path,in_col_name):
    """Wrapper to fix NANs in entrex IDs (IEEE-int has no NA)."""
    # here I have to check if there is a column named [in_col_name], and maybe NA check !
    df = pd.read_csv(tsv_file_path, sep="\t", dtype={in_col_name: str}, skipinitialspace=True)

    # remove gene version in ensembl
    df[in_col_name] = df[in_col_name].replace(to_replace=r'\.[0-9]+', value='', regex=True)

    return df

def convert_gene_ids(gene_ids, source, target, in_col_name, out_col_name, skip_missing, act_on_duplicated, tsv_conv_file_path):
    """Myege wrapper."""
    renamer = {
            "uniprot_id" : "uniprot.Swiss-Prot",
            "ensembl_id" : "ensembl.gene",
            "ensembl.transcript" : "ensembl.transcript",
            "entrez_id" : "entrezgene",
            "gene_symbol" : "symbol",
        }
    client = get_client("gene")
    genes_converted = client.querymany(
        gene_ids,
        scopes=renamer[source],
        fields=renamer[target],
        species="human",
        returnall=True,
        as_dataframe=True
    )
    ## return all is True: Extract duplicated and missing genes ##
    genes_dup = genes_converted["dup"]
    genes_missing = genes_converted["missing"]
    # Handle missing genes
    if genes_missing.size != 0:
        tsv_missing_file_path = tsv_conv_file_path.with_stem(tsv_conv_file_path.stem + "_missing")
        genes_missing.to_csv(tsv_missing_file_path, sep="\t", index=False)
        print(f"Wrote {tsv_missing_file_path} for missing genes.")
        if skip_missing==False :
            print(f"Note: {genes_missing.size} missing gene{'s' if genes_missing.size > 1 else ''}:{genes_missing.values.tolist()}")
            sys.stderr.write("Warning: Missing genes were found and they will be kept with an NA value.\n")
    # Handle duplicated genes
    if genes_dup.size != 0:
        print(f"Warning: {genes_dup.size} duplicated gene{'s' if genes_dup.size > 1 else ''}:{genes_dup.values.tolist()}")
        print(f"The action taken for these duplicated entries is:  {act_on_duplicated}")
        tsv_dup_file_path = tsv_conv_file_path.with_stem(tsv_conv_file_path.stem + "_dup")
        genes_dup.to_csv(tsv_dup_file_path, sep="\t", index=False)
        print(f"Wrote {tsv_dup_file_path} for duplicated genes")
    # keep missing values as NAN and reset index
    genes_converted = genes_converted["out"]
    if skip_missing==True :
        genes_converted = genes_converted.dropna(subset=renamer[target])
    genes_converted = genes_converted.reset_index(names=renamer[source])
    # act on duplicated genes here 
    if act_on_duplicated=="first_hit":
        genes_converted = genes_converted.drop_duplicates(subset=renamer[source])
    if act_on_duplicated=="remove_entry":
        # Remove duplicated rows from conv completely based on the in_col_name
        genes_converted = genes_converted.drop_duplicates(subset=renamer[source], keep=False)
    # Filter and Rename columns
    genes_converted = genes_converted.loc[:, [renamer[source], renamer[target]]]
    genes_converted.rename(columns={renamer[source]: in_col_name, renamer[target]: out_col_name}, inplace=True)
    return genes_converted

#https://click.palletsprojects.com/en/8.1.x/commands/ From Dani script
@click.group()
def main():
    """
        gene_id converter for dr-pipeline.
        This script allows you to convert gene identifiers from one format to another using the mygene library.
    """

    print(f"Running gene_id conversor {THIS_VERSION}")

@main.command()
@click.argument("tsv_file_path", type=click.Path(exists=True), required=True)
@click.option("--in_col_name", required=True, help="Name of the column containing the input gene identifiers.")
@click.option("--out_col_name", required=True, help="Name of the column to store the converted gene identifiers.")
@click.option("--in_id_format", type=click.Choice(["ensembl_id", "uniprot_id", "entrez_id", "gene_symbol"], case_sensitive=False), required=True, help="Format of the input gene identifiers.")
@click.option("--out_id_format", type=click.Choice(["entrez_id","gene_symbol"], case_sensitive=False), required=True, help="Desired format of the converted gene identifiers.")
@click.option('--skip_missing', default=False, help="Skip missing genes instead of exiting the script when encountered.")
@click.option("--act_on_duplicated", type=click.Choice(["first_hit","duplicate_entry", "remove_entry"], case_sensitive=False), default="duplicate_entry", help="Action to take when encountering duplicated gene identifiers.")
@click.argument("tsv_conv_file_path", type=click.Path(exists=False), required=True)
def translate(tsv_file_path, tsv_conv_file_path, in_col_name, out_col_name, in_id_format, out_id_format, skip_missing, act_on_duplicated):
    """
    Translate gene identifiers from one format to another.

    Args:

        tsv_file_path: Path to the TSV file containing gene identifiers.

        tsv_conv_file_path: Path to store the converted TSV file.

        in_col_name: Name of the column containing the input gene identifiers.

        out_col_name: Name of the column to store the converted gene identifiers.

        in_id_format: Format of the input gene identifiers.

        out_id_format: Desired format of the converted gene identifiers.

        skip_missing: Skip missing genes instead of keeping them with an NA value. (default: False)

        act_on_duplicated: Action to take when encountering duplicated gene identifiers (default: duplicate_entry).

    """
    print("Running mygene translation tool.")
    try:
        tsv_file_path = Path(tsv_file_path)
        tsv_conv_file_path = Path(tsv_conv_file_path)

        in_id_format = in_id_format.lower()
        out_id_format = out_id_format.lower()
        
        data = read_file_df(tsv_file_path,in_col_name)
        # Here check duplicated in the original data before proccessing!
        if any(data[in_col_name].duplicated()):
            sys.stderr.write("Note: The provided matrix contains duplicated genes.")
        ids = data[in_col_name].unique()
        genes_df = convert_gene_ids(ids, source=in_id_format, target=out_id_format, in_col_name=in_col_name, out_col_name=out_col_name, skip_missing=skip_missing, act_on_duplicated=act_on_duplicated,tsv_conv_file_path=tsv_conv_file_path)
        # Here add it to the original data with [out_col_name]
        conv_data_df = pd.merge(data, genes_df, on=in_col_name, how='left')
        if skip_missing==True :
            conv_data_df = conv_data_df.dropna(subset=out_col_name)
        conv_data_df.to_csv(tsv_conv_file_path, sep="\t", index=False)
        print(f"Wrote {tsv_conv_file_path}")
    except Exception as e:
        print(f"An error occurred: {str(e)}")
        sys.exit(1)  # Exit with non-zero status code to indicate failure    
@main.command()
@click.argument("tsv_file_path", type=click.Path(exists=True), required=True)
@click.option("--in_col_name", required=True, help="Name of the column containing the input gene identifiers.")
@click.option("--in_id_format", type=click.Choice(["ensembl_id", "uniprot_id", "entrez_id", "gene_symbol"], case_sensitive=False), required=True, help="Format of the input gene identifiers.")
@click.option('--returnall', default=False, help="Include all information for each gene, including duplicates and missing genes.")
@click.argument("tsv_conv_file_path", type=click.Path(exists=False), required=True)
def transtable(tsv_file_path, in_col_name, in_id_format, returnall, tsv_conv_file_path):
    """
    Generate a translation table for gene identifiers.

    Args:

        tsv_file_path: Path to the TSV file containing gene identifiers.

        in_col_name: Name of the column containing the input gene identifiers.

        in_id_format: Format of the input gene identifiers.

        returnall: Include all information for each gene, including duplicates and missing genes (default: False).

        tsv_conv_file_path: Path to store the generated translation table.

    """
    print("Translation table using mygene.")
    print(f"Returnall mode is {'on' if returnall else 'off'}")

    tsv_file_path = Path(tsv_file_path)
    tsv_conv_file_path = Path(tsv_conv_file_path)
    
    in_id_format = in_id_format.lower()

    data = read_file_df(tsv_file_path,in_col_name)
    gene_ids = data[in_col_name].unique()
    renamer = {
        "uniprot_id" :"uniprot.Swiss-Prot",
        "ensembl_id" : "ensembl.gene",
        "ensembl_lst" : "ensembl",
        "entrez_id" : "entrezgene",
        "symbol_id" : "symbol"
    }
    client = get_client("gene")
    genes_converted = client.querymany(
        gene_ids,
        scopes=renamer[in_id_format],
        fields=list(renamer.values()),
        species="human",
        returnall=returnall,
        as_dataframe=True,
    )
    if returnall:
        genes_dup = genes_converted["dup"]
        genes_missing = genes_converted["missing"]
        genes_converted = genes_converted["out"]
    cols_query = genes_converted.columns.isin(list(renamer.values()))
    genes_converted = genes_converted.loc[:, cols_query]
    genes_converted = genes_converted.rename(columns={v: k for k, v in renamer.items()})
    # Here add a validator
    print(genes_converted)
    
    # Here add it to the original data with [out_col_name]
    genes_converted.to_csv(tsv_conv_file_path, sep="\t", index=True)
    print(f"Wrote {tsv_conv_file_path}")
    if returnall:
        print(f"Warning: There are {genes_dup.size} duplicated genes and {genes_missing.size} missing genes.")
        if genes_dup.size != 0:
            tsv_dup_file_path = tsv_conv_file_path.with_stem(tsv_conv_file_path.stem + "_dup")
            genes_dup.to_csv(tsv_dup_file_path, sep="\t", index=False)
            print(f"Wrote {tsv_dup_file_path} for duplicated genes")
        if genes_missing.size != 0:
            tsv_missing_file_path =  tsv_conv_file_path.with_stem(tsv_conv_file_path.stem +"_missing")
            genes_missing.to_csv(tsv_missing_file_path, sep="\t", index=False)
            print(f"Wrote {tsv_missing_file_path} for missing genes.")
if __name__ == "__main__":
    main()
