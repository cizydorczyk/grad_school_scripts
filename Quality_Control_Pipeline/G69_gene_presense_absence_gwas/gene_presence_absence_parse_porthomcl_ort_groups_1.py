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
import itertools

script, isolate_list, orthologs_file, outputdir = argv

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

# Parse input files:
isolates_list = []
with open(isolate_list, 'r') as infile1:
    for line in infile1:
        isolates_list.append(line.strip())

gene_families = []
with open(orthologs_file, 'r') as infile2:
    for line in infile2:
        gene_families.append(line.strip().split('\t'))

# Make output directory:
outpath = outputdir + "fishers_exact_test_gene_presence_absence"
os.makedirs(outpath)

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

# Convert binary data to dataframe:
df1 = pd.DataFrame.from_dict(gene_presence_absence_dict, orient="index")
# Name dataframe columns (by family number - same order as in orthologs file)""
fam_names = range(1, len(gene_families) + 1)
df1.columns = fam_names
# Add phenotypes column to dataframe:
df1['Phenotype'] = pd.Series(eradication_dict)




# Run fishers exact test:
to_write_all_contingency_tables = []
to_write_good_contingency_tables = []
to_correct_contingency_tables = set()

# Get raw pvalues using Fisher's exact test:
rawpvalues = []
rawoddsratios = []
good_fam_names = []

# To_drop not currently used for anything, but it's a list of fam names that are not testable:
to_drop =  []

# Identify testable gene families:
for i in fam_names:
    tab = pd.crosstab(df1["Phenotype"] == '1', df1[i] == '1')
    to_write_all_contingency_tables.append(tab.to_string())

    if tab.shape == (2,2):
        temp1 = tab.values.tolist()
        temp2 = temp1[0] + temp1[1]
        if temp2.count(0) == 1:
            to_drop.append(i)
        else:
            # Run fisher's exact test:
            oddsratio, p_value = fisher_exact(tab)
            rawpvalues.append(p_value)
            rawoddsratios.append(oddsratio)
            good_fam_names.append(i)

            to_write_good_contingency_tables.append(tab.to_string())

            # Get effective number of contingency tables (for p-value correction):
            temp3 = map(str, temp2)
            temp4 =  '.'.join(temp3)
            to_correct_contingency_tables.add(temp4)
    else:
        to_drop.append(i)



df2 = df1.drop(to_drop, axis=1)

# Get number of gene profiles:
gene_profile_dict = df2.to_dict(orient="list")

set_gene_profiles = set()
for i in gene_profile_dict:
    set_gene_profiles.add(''.join(gene_profile_dict[i]))

# Get pvalue correction values and print to terminal:
strict_bonferroni_correction = len(good_fam_names)
gene_fam_correction = len(set_gene_profiles)
contingency_table_correction = len(to_correct_contingency_tables)

print "Number of total gene families: ", len(gene_families)
print "Number of gene families: ", strict_bonferroni_correction
print "Number of gene profiles: ", gene_fam_correction
print "Number of contingency table profiles: ", contingency_table_correction

# Correct pvalues:
strict_bonferroni_pvalues = []
gene_fam_pvalues = []
contingency_table_pvalues = []

for i in rawpvalues:
    strict_bonferroni_pvalues.append(i*strict_bonferroni_correction)
    gene_fam_pvalues.append(i*gene_fam_correction)
    contingency_table_pvalues.append(i*contingency_table_correction)

# Write output:
to_write = []
to_write.append('Family_Number' + '\t' + 'Raw_Odds_Ratio' + '\t' + 'Raw_p-value' + '\t' + 'Strict_BF_p-value' + '\t' + 'Gene_Profile_BF_p-value' + '\t' + 'Contingency_Table_BF_p-value')
for i, j, k, l, m, n in itertools.izip(good_fam_names, rawoddsratios, rawpvalues, strict_bonferroni_pvalues, gene_fam_pvalues, contingency_table_pvalues):
    to_write.append(str(i) + '\t' + str(j) + '\t' + str(k) + '\t' + str(l) + '\t' + str(m) + '\t' + str(n))

with open(os.path.join(outpath, "fishers_exact_results.txt"), 'w') as outfile1:
    outfile1.write('\n'.join(to_write))

with open(os.path.join(outpath, "good_contingency_tables.txt"), 'w') as outfile2:
    outfile2.write('\n'.join(to_write_good_contingency_tables))

with open(os.path.join(outpath, "raw_contingency_tables.txt"), 'w') as outfile2:
    outfile2.write('\n'.join(to_write_all_contingency_tables))

df1_output = os.path.join(outpath, "raw_gene_presence_absence_matrix.txt")
df2_output = os.path.join(outpath, "filtered_gene_presence_absence_matrix.txt")
df1.to_csv(path_or_buf=df1_output, sep='\t')
df2.to_csv(path_or_buf=df2_output, sep='\t')

    # output_contingency.append(tab.to_string())
    # # Get number of identical tables:
    # temp1 = tab.values.tolist()
    # temp2 = temp1[0]
    # temp3 = temp1[1]
    # temp4 = temp2 + temp3
    # # Convert all ints in temp4 list to string:
    # temp5 = map(str, temp4)
    # temp6 = ''.join(temp5)
    # output_contingency_strs.add(temp6)




    # Run fisher's exact test:
    # try:
    #     oddsratio, p_value = fisher_exact(tab)
    #     rawpvalues.append(p_value)
    #     rawoddsratios.append(oddsratio)
    # except ValueError:
    #     oddsratio = "NA"
    #     p_value = "NA"
    #     rawpvalues.append(p_value)
    #     rawoddsratios.append(oddsratio)





# # Get pvalue corrections:
# # print len(fam_names)
# strict_bonferroni_correction = len(fam_names)
# # print len(set_gene_profiles)
# gene_fam_correction = len(set_gene_profiles)
# # print len(output_contingency_strs)
# contingency_tables_correction = len(output_contingency_strs)
#
# # Correct pvalues:
# strict_bonferroni_pvalues = []
# for i, j in itertools.izip(rawpvalues, rawoddsratios):
#     if i == "NA":
#         strict_bonferroni_pvalues.append("NA")
#     elif j == 'inf':
#         strict_bonferroni_pvalues.append("NA")
#     elif j == '0.0':
#         strict_bonferroni_pvaules.append("NA")
#     else:
#         corrected_pvalue = i*strict_bonferroni_correction
#         strict_bonferroni_pvalues.append(corrected_pvalue)
#
# gene_fam_pvalues =[]
# for i, j in itertools.izip(rawpvalues, rawoddsratios):
#     if i == "NA":
#         gene_fam_pvalues.append("NA")
#     elif j == 'inf':
#         gene_fam_pvalues.append("NA")
#     elif j == '0.0':
#         gene_fam_pvalues.append("NA")
#     else:
#         corrected_pvalue = i*strict_bonferroni_correction
#         gene_fam_pvalues.append(corrected_pvalue)
#
# contingency_tables_pvalues = []
# for i, j in itertools.izip(rawpvalues, rawoddsratios):
#     if i == "NA":
#         contingency_tables_pvalues.append("NA")
#     elif j == 'inf':
#         contingency_tables_pvalues.append("NA")
#     elif j == '0.0':
#         contingency_tables_pvalues.append("NA")
#     else:
#         corrected_pvalue = i*strict_bonferroni_correction
#         contingency_tables_pvalues.append(corrected_pvalue)
#
# # strict_bonferroni_pvalues = []
# # for i, j in itertools.izip(rawpvalues, rawoddsratios):
# #     if i != "NA" and (j != "inf" or j != '0.0'):
# #         corrected_pvalue = i*strict_bonferroni_correction
# #         strict_bonferroni_pvalues.append(corrected_pvalue)
# #     else:
# #         strict_bonferroni_pvalues.append("NA")
# #
# # gene_fam_pvalues = []
# # for i, j in itertools.izip(rawpvalues, rawoddsratios):
# #     if i != "NA" and (j != "inf" or j != '0.0'):
# #         corrected_pvalue = i*gene_fam_correction
# #         gene_fam_pvalues.append(corrected_pvalue)
# #     else:
# #         gene_fam_pvalues.append("NA")
# #
# # contingency_tables_pvalues = []
# # for i, j in itertools.izip(rawpvalues, rawoddsratios):
# #     if i != "NA" and (j != "inf" or j != '0.0'):
# #         corrected_pvalue = i*contingency_tables_correction
# #         contingency_tables_pvalues.append(corrected_pvalue)
# #     else:
# #         contingency_tables_pvalues.append("NA")
#

#
