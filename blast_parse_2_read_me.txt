This script takes both the blast output (-outfmt "7 qseqid sseqid pident length
evalue sscinames") and CLC Genomics Workbench assembly statistics file (export
contig name, contig length, coverage, etc.) as a tsv (tab delimited) file. It
parses these files and creates a new tsv file with columns corresponding to
contig name, length, coverage, inner color, outer color, and what it blasts to.

The point of this new tsv file is that it can be used in R to graph contig
length vs. coverage, making each contig (dot) that BLASTs to Pa blue and each
one that doesn't yellow.

It takes three arguments:

    inputblastfile: the input blast file (formatted as above)

    inputcontigstats: CLC genomics stats tsv Filename

    output_tsv: output tsv file read for R
