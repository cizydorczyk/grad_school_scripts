This script takes 4 arguments:

    snp file created using 3_aligner_snp_call_v1.py
    TAB genome annotation from pseudomonas.com
    GTF annotation file with intergenic regions from pseudomonas.com
    output annotated snp file

The script basically takes a snp file and for each snp, identifies whether it
falls in a gene or intergenic region, and returns info on that snp based on its
category.

Intergenic snps are classified in one of three categories: regulatory region left
(RRL), regulatory region none (RRN), or regulatory region right (RRR).

It basically parses each input annotation file (except output...duh) and creates
dictionaries of either class objects or lines, each indexed by the ranges of
each gene/intergenic region.

It parses the input snps file also as a dictionary, with each key being the snp
position.

If the snp position falls within a gene range in the genes dictionary, info is
pulled from that dictionary for the snp. Else if the snp position falls in an
intergenic region, info is pulled from the integenic regions dictionary for the
snp. It is assummed all snps fall into either one of these dictionaries or the
other. Aka if a snp doesn't fall in either, something is wrong with the input
file, and the script will spit up some error (by default; I didn't add error
reporting capability).

See script for comments.

#### MAKE SURE THERE ARE NO BLANK LINES AT THE END OF THE FILE ####
