## Everything is run by Snakemake

## Mixes
`Rscript scatterplot.R [Miseq|454] [GEO|20-20|60-10|80-5|90_2.5]` creates tables and pdf figures.

## Patients

Script within Snakefile is run, then table for patient 2 is obtained with

    tr ',' '-' < patient_2_merged.tsv | tr -d ' ' | tr '\t' ','

The output is then copied into
[Tables Generator](http://www.tablesgenerator.com/markdown_tables)
