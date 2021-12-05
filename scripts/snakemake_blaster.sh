#!/bin/bash

## Prerequisites: install tools:
# conda create --name taiga -c bioconda  blast biopython snakemake
## Download BLAST nr database from https://ftp.ncbi.nlm.nih.gov/blast/db/
## Extract taxids of interest, e.g., for bacteria (the script is available from the BLAST installation package):
# get_species_taxids.sh -t 2 > bacterial.ids
## Adjust snakemake file with proper path to database, taxids, input and output directories
# vim snakemake_blaster.snakemakefile

THREADS=10

snakemake --snakefile snakemake_blaster.snakemakefile -k -j $THREADS --rerun-incomplete
