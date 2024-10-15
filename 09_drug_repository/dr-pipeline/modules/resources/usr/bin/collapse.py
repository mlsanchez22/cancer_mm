#!/usr/bin/env python
# -*- coding: utf-8 -*-
"""
Collapse RNAseq counts by column unique IDs.
"""

import click
import pandas as pd

THIS_VERSION = 1.0


@click.command()
@click.argument(
    'input_file',
    required=True, 
    type=click.Path(exists=True))
@click.argument(
    'output_file',
    required=False, 
    type=click.Path(exists=False), 
    default="output.txt")
@click.option(
    '--how', 
    type=click.Choice(['max', 'mean', 'sum'], case_sensitive=False),
    default='sum',
    show_default=True,
    help='Collapse the lines with duplicate identifiers by applying the sum, the mean, or the maximum value.')
@click.option(
    "--id_column",
    required=True,
    help="Column name used as unique identifiers.")
@click.option(
    "--skiprows",
    default=0,
    help="Skip n rows in the input file.")
def main(input_file, output_file, how, id_column, skiprows):
    """RNAseq collapser."""

    print(f"Running RNAseq collapser {THIS_VERSION}")

    print(f"Loading input file {input_file}")
    df = pd.read_csv(
        input_file,
        sep="\t",
        skiprows=skiprows,
        index_col=[0,1]
    )
    print(f"...the input data set contains {df.shape[0]} rows and {df.shape[1]} columns")

    print(f"Grouping rows by \"{id_column}\" column")
    g = df.groupby(by=id_column)
    filtered = g.filter(lambda x: len(x) >= 2)
    print(f"...the collapsed dataset contains {g.ngroups} rows.")
    #print(f"...a total of {filtered[id_column].unique().size} IDs are duplicated in the original dataset.")

    print(f"Writting {how} collapsed results to {output_file}")
    g.agg(how).to_csv(output_file, sep="\t", index=True)

    print(f"Writting filtered rows file")
    filtered.to_csv("duplicated.txt", sep="\t", index=True)

    print("Done!")


if __name__ == "__main__":
    main()
