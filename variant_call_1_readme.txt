This script takes vcf files, filters them, and identifies high quality variant
calls (SNPs).

This script takes argumnents:

    Log file (text file to contain info from running the script)
    List (text file, one full file path per line) of vcf files
    Output fasta multi-sequence alignment file

This script follows the logical approach outlined in variant_call_1_pseudocode.txt.

Below is a walkthrough of the script.

The script calls the following three functions on every file in the input list:

    parse_vcf
    indel_position_list
    hq_positions_list

parse_vcf() takes a list of vcf files (full paths) and parses them, creating a
VcfVariant class objects (class defined in script) for each variant.

indel_position_list then goes through the vcf variant objects list and notes the
position of variants that are indels (have a reference or alternate base > 1bp).

hq_positions_list uses the output from the previous two functions to identify
high quality variants and notes their positions, as well as the reference base
at that position.

    It filters variants by the following criteria:

        DP (depth) > 20
        Quality > 30
        Not within 150 bp of the ends of the reference genome or an indel (
            based on  indel positions list from function 2)
        Alternate reads in the f and r directions must number at least 3
        Minimum % difference between alternate f and r reads is 20%

    This function doesn't return anything, but instead modifies global lists/
    dictionaries (hq_positions_list1 and hqrefdict) defined prior to the
    function call.

The script then parses the vcf files and identifies indel positions (same as above)
and then calls the "list_of_SNPs()" function, which, using the positions from
hq_positions_list() function, relaxes filtering criteria and creates a list of
vcf variant class objects at these positions.

    Relaxed criteria means:

        DP > 15
        Quality > 20

The script then calls whether a variant is a SNP (>80% of reads support alternate call),
or a reference call (< 20% of reads support alternate base call). If 20-80% of
reads support alternate call, the position is called as ambiguous ("N"). If a
high quality position does not have a variant in a particular file, the reference
base is called at that position.

The script runs the above functions for every file in the input list, after which
it creates a sequence for each file consisting of base calls (reference or alternate)
for each position in the high quality positions list, and outputs it in fasta format
in a single file for all vcf files.

    NOTE that it creates the sequences in order of positions (ie. lowest positions
    first).

    This output fasta file is considered a multiple sequence alignment, and is
    ready for downstream analysis (eg. phylogenetic).
