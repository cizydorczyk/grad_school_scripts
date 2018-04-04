from sys import argv
import os
import itertools

script, raw_fisher_results, orthologs_file = argv

class FisherOutput(object):
    def __init__(self, genefam, oddsratio, pvalue, strict_correction, gene_fam_correction,  cont_table_correction, record):
        self.genefam = genefam
        self.oddsratio = oddsratio
        self.pvalue = pvalue
        self.strict_correction = strict_correction
        self.gene_fam_correction = gene_fam_correction
        self.cont_table_correction = cont_table_correction
        self.record = record

fisher_objects = []
with open(raw_fisher_results, 'r') as infile1:
    for line in infile1:
        if not line.startswith("Family_Number"):
            temp_line = line.strip().split('\t')
            fisher_objects.append(FisherOutput(int(temp_line[0]), float(temp_line[1]), float(temp_line[2]), float(temp_line[3]), float(temp_line[4]), float(temp_line[5]), '\t'.join(temp_line)))

significant_cont_table_objects = []
significant_gene_fam_objects = []
significant_strict_objects = []


for i in fisher_objects:
    if i.cont_table_correction < 0.05:
        significant_cont_table_objects.append(i)
    if i.gene_fam_correction < 0.05:
        significant_gene_fam_objects.append(i)
    if i.strict_correction < 0.05:
        significant_strict_objects.append(i)

print len(significant_strict_objects)
print len(significant_gene_fam_objects)
print len(significant_cont_table_objects)

sig_cont_table_records = [i.record for i in significant_cont_table_objects]
sig_genefam_records = [i.record for i in significant_gene_fam_objects]
sig_strict_records = [i.record for i in significant_strict_objects]

cwd = os.getcwd()
with open(os.path.join(cwd, "cont_table__profile_corr_sig_results.txt"), 'w') as outfile3:
    outfile3.write('Family_Number' + '\t' + 'Raw_Odds_Ratio' + '\t' + 'Raw_p-value' + '\t' + 'Strict_BF_p-value' + '\t' + 'Gene_Profile_BF_p-value' + '\t' + 'Contingency_Table_BF_p-value' + '\n')
    outfile3.write('\n'.join(sig_cont_table_records))

with open(os.path.join(cwd, "gene_profile_corr_sig_results.txt"), 'w') as outfile2:
    outfile2.write('Family_Number' + '\t' + 'Raw_Odds_Ratio' + '\t' + 'Raw_p-value' + '\t' + 'Strict_BF_p-value' + '\t' + 'Gene_Profile_BF_p-value' + '\t' + 'Contingency_Table_BF_p-value' + '\n')
    outfile2.write('\n'.join(sig_genefam_records))

with open(os.path.join(cwd, "total_num_fam_corr_sig_results.txt"), 'w') as outfile1:
    outfile1.write('Family_Number' + '\t' + 'Raw_Odds_Ratio' + '\t' + 'Raw_p-value' + '\t' + 'Strict_BF_p-value' + '\t' + 'Gene_Profile_BF_p-value' + '\t' + 'Contingency_Table_BF_p-value' + '\n')
    outfile1.write('\n'.join(sig_strict_records))

with open(os.path.join(cwd, "summary.txt"), 'w') as outfile4:
    outfile4.write("Number of significant results using BF correction for total number of gene families: " + str(len(significant_strict_objects)) + '\n')
    outfile4.write("Number of significant results using BF correction for number of gene profiles: " + str(len(significant_gene_fam_objects)) + '\n')
    outfile4.write("Number of significant results using BF correction for number of contingency table profiles: " + str(len(significant_cont_table_objects)))

################
# Get gene families that are significant at each level and write to file:

gene_families = []
with open(orthologs_file, 'r') as infile2:
    for line in infile2:
        gene_families.append(line.strip().split('\t'))

gene_families_dict = {}
for i, j in itertools.izip(range(1, len(gene_families) + 1), gene_families):
    gene_families_dict[i] = j

strict_fams = []
genefam_fams = []
cont_table_fams = []

for i in gene_families_dict:
    if i in [j.genefam for j in significant_strict_objects]:
        strict_fams.append('\t'.join(gene_families_dict[i]))
    if i in [k.genefam for k in significant_gene_fam_objects]:
        genefam_fams.append('\t'.join(gene_families_dict[i]))
    if i in [l.genefam for l in significant_cont_table_objects]:
        cont_table_fams.append('\t'.join(gene_families_dict[i]))

with open(os.path.join(cwd, "total_num_fam_corr_gene_families.txt"), 'w') as outfile5:
    outfile5.write('\n'.join(strict_fams))

with open(os.path.join(cwd, "gene_profile_corr_gene_familes.txt"), 'w') as outfile6:
    outfile6.write('\n'.join(genefam_fams))

with open(os.path.join(cwd, "cont_table_corr_gene_families"), 'w') as outfile7:
    outfile7.write('\n'.join(cont_table_fams))
