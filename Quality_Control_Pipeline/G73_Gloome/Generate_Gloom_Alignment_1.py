from sys import argv
import os

script, isolate_list, input_families = argv

# Parse input files:
isolates_list = []
with open(isolate_list, 'r') as infile1:
    for line in infile1:
        isolates_list.append(line.strip())

gene_families = []
with open(input_families, 'r') as infile2:
    for line in infile2:
        gene_families.append(line.strip().split('\t'))

# Get working directory:
cwd = os.getcwd()

# Initialize gene presence/absence dictionary, with a key for each isolate:
gene_presence_absence_dict = {}
for i in isolates_list:
    gene_presence_absence_dict[i] = []

# Convert gene presence/absence to binary data:
for genefam in gene_families:
    # Get isolates present in family:
    isolates_present = set()

    for i in genefam:
        isolate_number = i.strip().split('|')[0].split('_')[0]
        isolates_present.add(isolate_number)

    # Record presence/absence of each gene family for each isolate:
    for isolate in isolates_list:
        if isolate in isolates_present:
            gene_presence_absence_dict[isolate].append('1')
        else:
            gene_presence_absence_dict[isolate].append('0')

to_write = []
for key in gene_presence_absence_dict:
	to_write.append(">" + str(key) + '\n' + ''.join(gene_presence_absence_dict[key]))

# Write gene presence/absence alignment to file:
with open(os.path.join(cwd, "gene_fam_corr_pres_abs_aln.txt"), 'w') as outfile1:
	outfile1.write('\n'.join(to_write))
