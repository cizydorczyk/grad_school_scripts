'''I also wrote a python script, "parse_porthomcl_ortho_group_output_1.py", to parse the ortholog family output from Porthomcl.
It takes the "8.all.ort.group" output file from PorthoMCL, an isolate list file, and a text file with the complete paths to all
isolate .faa files (from Prokka). It outputs a couple summary files (number of core and accessory families, and counts of how many
isolates each family is found in (# of isolates : # of families present in that many isolates). It also outputs two folders: one
containing fasta files for each core family (one file per family), and one for the accessory families (one file per family).
NOTE that these are NOT alignments, just the fasta sequences that make up each family. There may be (and are!) multiple genes in
each family from the same isolate - this script just clumps them all together.'''

from sys import argv
from Bio import SeqIO
import itertools
import os

script, isolate_list, orthologs_file,  genome_faa_files_list, outputdir = argv

isolates_list = []
with open(isolate_list, 'r') as infile1:
    for line in infile1:
        isolates_list.append(line.strip())

gene_families = []
with open(orthologs_file, 'r') as infile2:
    for line in infile2:
        gene_families.append(line.strip().split('\t'))

core_families = []
accessory_families = []

family_sizes = {}
for i in range(1, 156):
    family_sizes[i] = 0

missing_isolates = {}
for i in range(1, 156):
    missing_isolates[i] = []

outpath = outputdir + "core_fam_missing_counts"
os.makedirs(outpath)

# loop goes here...
for number in range(1,156):
    for genefam in gene_families:
        isolates_present = set()
        for i in genefam:
            isolate_number = i.strip().split('|')[0].split('_')[0]
            isolates_present.add(isolate_number)

        if len(isolates_present) == number:
            for i in isolates_list:
                if i not in isolates_present:
                    missing_isolates[number].append(i)

        family_sizes[len(isolates_present)] += 1

    missing_counts = {}
    set154 = set(missing_isolates[number])
    for i in set154:
        missing_counts[i] = 0

    for j in missing_isolates[number]:
        missing_counts[j] += 1

    print missing_counts

    with open(os.path.join(outpath, str(number) + "_counts.txt"), 'w') as outfile1:
        to_write = []
        for i in sorted(missing_counts):
            to_write.append(str(i) + '\t' + str(missing_counts[i]))
        outfile1.write('\n'.join(to_write))
