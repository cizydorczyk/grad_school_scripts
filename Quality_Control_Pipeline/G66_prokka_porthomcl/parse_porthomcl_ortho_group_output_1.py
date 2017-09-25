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

for genefam in gene_families:
    isolates_present = set()
    for i in genefam:
        isolate_number = i.strip().split('|')[0].split('_')[0]
        isolates_present.add(isolate_number)
    if len(isolates_present) == 155:
        core_families.append(genefam)
    elif len(isolates_present) != 155:
        accessory_families.append(genefam)

    family_sizes[len(isolates_present)] += 1

print len(core_families)
print len(accessory_families)
print family_sizes

with open(genome_faa_files_list, 'r') as infile3:
    file_list = []
    for line in infile3:
        file_list.append(line.strip())

isolate_gene_dicts =  {}
for j,i in itertools.izip(isolates_list, file_list):
    isolate_id = str(j) + "_pae"
    isolate_gene_dicts[isolate_id] = SeqIO.to_dict(SeqIO.parse(i.strip(), "fasta"))

# Collect core genes for each core family:
core_fam_number = range(1, len(core_families) + 1)
core_fam_fasta_records = {}

for corenumber, genefam in itertools.izip(core_fam_number, core_families):
    genefam_records = []
    for gene in genefam:
        geneisolate = gene.strip().split('|')[0]
        geneid = gene.strip().split('|')[1]

        try:
            gene_fasta_entry = isolate_gene_dicts[geneisolate][geneid]
            genefam_records.append(gene_fasta_entry)

        except KeyError:
            geneisolate2 =  list(geneisolate)
            geneisolate3 =  geneisolate2[0] + geneisolate2[1] + geneisolate[2] + '-' + geneisolate2[3] + "_pae"
            gene_fasta_entry = isolate_gene_dicts[geneisolate3][geneid]
            genefam_records.append(gene_fasta_entry)

    core_fam_fasta_records[corenumber] = genefam_records
print len(core_fam_fasta_records)

# Collect accessory genes for each accessory gene family:
acc_fam_number = range(1, len(accessory_families) + 1)
acc_fam_fasta_records = {}

for accnumber, genefam in itertools.izip(acc_fam_number, accessory_families):
    genefam_records = []
    for gene in genefam:
        geneisolate = gene.strip().split('|')[0]
        geneid = gene.strip().split('|')[1]

        try:
            gene_fasta_entry = isolate_gene_dicts[geneisolate][geneid]
            genefam_records.append(gene_fasta_entry)

        except KeyError:
            geneisolate2 =  list(geneisolate)
            geneisolate3 =  geneisolate2[0] + geneisolate2[1] + geneisolate[2] + '-' + geneisolate2[3] + "_pae"
            gene_fasta_entry = isolate_gene_dicts[geneisolate3][geneid]
            genefam_records.append(gene_fasta_entry)

    acc_fam_fasta_records[accnumber] = genefam_records
print len(acc_fam_fasta_records)

############# Write output files ##############

outpath = outputdir + "core_fam_faa_records"
os.makedirs(outpath)

outpath2 = outputdir + "acc_fam_faa_records"
os.makedirs(outpath2)

for corenumber, genefam in itertools.izip(core_fam_number, core_fam_fasta_records):
    output_file = "core_fam_" + str(corenumber) + ".faa"
    SeqIO.write(core_fam_fasta_records[genefam], os.path.join(outpath, output_file), "fasta")

for accnumber, genefam in itertools.izip(acc_fam_number, acc_fam_fasta_records):
    output_file = "acc_fam_" + str(accnumber) + ".faa"
    SeqIO.write(acc_fam_fasta_records[genefam], os.path.join(outpath2, output_file), "fasta")

with open(os.path.join(outputdir, "summary.txt"), 'w') as outfile2:
    outfile2.write("Number of core gene families: " + str(len(core_families)) + '\n' +
    "Number of accessory gene families: " + str(len(accessory_families)))

with open(os.path.join(outputdir, "isolates_per_fam_summary.txt"), 'w') as outfile3:
    to_write = []
    for i in sorted(family_sizes):
        to_write.append(str(i) + '\t' + str(family_sizes[i]))
    outfile3.write("Number_Isolates" + '\t' + 'Number_of_Families' + '\n')
    outfile3.write('\n'.join(to_write))
