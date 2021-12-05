# About

This repository include scripts for analysis of biosynthetic gene clusters (BGCs) in large metagenome assemblies. All scripts were created within the [Chernevaya Taiga project](https://cab.spbu.ru/grant/rsf-19-16-00049/) funded by the Russian Science Foundation (grant 19-16-00049). The repository also contains results of the scripts application to four hybrid (illumina + ONT) assemblies of various Chernevaya Taiga soil samples.

# Suggested analysis pipeline
1. Raw sequencing data quality control and quality trimming if needed. Recommended software: FastQC, Trimmomatic.
2. Metagenome assembly. Recommended software: metaSPAdes (for illumina and illumina+ONT), metaFlye (for ONT only data)
3. BGC identification. Recommended software: antiSMASH.
4. BGC summary statistics. Script: `scripts/antismash_json_parser.py`
5. BGC core and additional biosynthetic gene sequence extraction. Script: `scripts/antismash_json_parser.py`
6. BLASTing of the extracted genes against NCBI databases. Script: `scripts/snakemake_blaster.sh`
7. Summarizing BLAST output. Script: `scripts/blast_xml_parser.py`

# Contact
Please send your questions, feedback, bug reports to Alexey Gurevich: <aleksey.gurevich@spbu.ru> or post a GitHub issue about them.