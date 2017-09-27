'''I also wrote a python script, "parse_porthomcl_ortho_group_output_1.py", to parse the ortholog family output from Porthomcl.
It takes the "8.all.ort.group" output file from PorthoMCL, an isolate list file, and a text file with the complete paths to all
isolate .faa files (from Prokka). It outputs a couple summary files (number of core and accessory families, and counts of how many
isolates each family is found in (# of isolates : # of families present in that many isolates). It also outputs two folders: one
containing fasta files for each core family (one file per family), and one for the accessory families (one file per family).
NOTE that these are NOT alignments, just the fasta sequences that make up each family. There may be (and are!) multiple genes in
each family from the same isolate - this script just clumps them all together.'''


# Note that gene families are in the same order as they appear in the orthologs file.
# This means that the first family in the dataframe is equal to the first family in
# the orthologs file.

from sys import argv
import os
import numpy as np
import pandas as pd
from scipy.stats import fisher_exact

script, isolate_list, orthologs_file,  genome_faa_files_list, outputdir = argv

eradication_dict = {'6':'0', '7':'0', '8':'0', '10':'0', '14':'0', '17':'0', '19':'0',
'21':'0', '23':'0', '30':'0', '36':'0', '37':'0', '53':'0', '97':'0', '130':'0', '152':'0', '160':'0',
'180':'0', '196':'0', '200':'0', '217':'0', '225':'0', '235':'0', '251':'0', '257':'1',
'258':'1', '259':'0', '263':'0', '265':'0', '270':'0', '273':'0', '276':'0', '282':'0',
'275':'0', '278':'0', '280':'0', '288':'0', '290':'0', '293':'0', '298C':'1', '303':'0',
'306':'0', '310':'0', '315':'0', '317':'0', '318':'0', '323':'0', '325':'0', '326':'0',
'327':'0', '330':'1', '331':'1', '332':'1', '335':'1', '336':'1', '337':'1', '340':'1',
'341':'1', '342':'1', '358':'0', '359':'0', '360':'0', '366':'0', '367':'0', '369':'0',
'375':'1', '379':'1', '380':'1', '385':'1', '388':'1', '390':'0', '392':'0', '393':'0',
'395':'0', '400':'0', '401':'0', '402':'0', '404C':'0', '405C':'0', '406C':'0', '409':'0',
'410C':'0', '411C':'0', '412':'0', '417C':'0', '418C':'0', '419C':'0', '420N2':'0',
'427':'0', '428':'0', '429':'0', '430':'0', '431':'0', '432':'0', '433':'0', '446':'0',
'448':'0', '449':'0', '450':'0', '453':'0', '455':'0', '457':'0', '459C':'0', '465':'1',
'471':'0', '472':'0', '473':'0', '475':'0', '476':'0', '479':'0', '480':'0', '487':'0',
'501':'0', '503':'1', '504':'1', '505':'1', '506':'1', '507':'0', '508':'0', '509':'0',
'510':'0', '511':'0', '512':'0', '513':'0', '514':'0', '521':'0', '525':'0', '526':'0',
'527':'0', '539':'0', '540':'0', '5412':'0', '544':'0', '551':'1', '552':'1', '553C':'1',
'557':'0', '5592':'0', '562':'0', '563':'0', '566':'1', '567':'1', '568':'1', '569':'0',
'570':'0', '571':'0', '5721':'1', '5732':'1', '5733':'1', '575':'0', '577':'0', '578':'0',
'580':'1', '581':'1', '583':'1'}

isolates_list = []
with open(isolate_list, 'r') as infile1:
    for line in infile1:
        isolates_list.append(line.strip())

gene_families = []
with open(orthologs_file, 'r') as infile2:
    for line in infile2:
        gene_families.append(line.strip().split('\t'))

outpath = outputdir + "fishers_exact_test_gene_presence_absence"
os.makedirs(outpath)

# Initialize gene presence/absence dictionary, with a key for each isolate:
gene_presence_absence_dict = {}
for i in isolates_list:
    gene_presence_absence_dict[i] = []

# Main algorithm:
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

df1 = pd.DataFrame.from_dict(gene_presence_absence_dict, orient="index")
fam_names = range(1, len(gene_families) + 1)
df1.columns = fam_names
df1['Phenotype'] = pd.Series(eradication_dict)

# Run fishers exact test:
output_list = []
output_contingency = []

for i in fam_names:
	tab = pd.crosstab(df1["Phenotype"] == '1', df1[i] == '1')
	output_contingency.append(tab.to_string())
	print tab
	try:
		oddsratio, p_value = fisher_exact(tab)
		output_list.append(str(i) + '\t' + str(oddsratio) + '\t' + str(p_value))
	except ValueError:
		oddsratio = "NA"
		p_value = "NA"
		output_list.append(str(i) + '\t' + oddsratio + '\t' + p_value)

with open(os.path.join(outpath, "fishers_exact_naive_results.txt"), 'w') as outfile1:
	outfile1.write('SNP' + '\t' + 'oddsratio' + '\t' + 'p_value' + '\n')
	outfile1.write('\n'.join(output_list))

with open(os.path.join(outpath, "fishers_exact_contingency_tables.txt"), 'w') as outfile2:
	outfile2.write('\n'.join(output_contingency))
