from sys import argv
import numpy as np
import pandas as pd
import itertools

script, inputtxt, outputfile1 = argv

df = pd.read_table(inputtxt)
df = df.transpose()
#df.drop("Total_Count", axis=0, inplace=True)
df.columns = df.iloc[0]
df.drop("Gene", axis=0, inplace=True)

df.to_csv(path_or_buf="/home/conrad/Data/h_notebook_files/H6_treewas_maf005_genes_uncollapsed/gene_presence_absence_maf005_transposed.txt", sep='\t')

df_dict = df.to_dict(orient="list")
df_list = list(df_dict)

guide_dict = {}
unique_profiles_list = []

combined_snp_positions = {}

for key, value in df_dict.iteritems():
    if value not in unique_profiles_list:
        unique_profiles_list.append(value)

for key, value in df_dict.iteritems():

    if unique_profiles_list.index(value) not in guide_dict:
        guide_dict[unique_profiles_list.index(value)] = [key]
    elif unique_profiles_list.index(value) in guide_dict:
        guide_dict[unique_profiles_list.index(value)].append(key)

    if unique_profiles_list.index(value) not in combined_snp_positions:
        combined_snp_positions[unique_profiles_list.index(value)] = value
    elif unique_profiles_list.index(value) in combined_snp_positions:
        continue

print combined_snp_positions[5]

with open(outputfile1, 'w') as outfile1:
    for i, j in itertools.izip(range(0, len(unique_profiles_list)), unique_profiles_list):
        outfile1.write(str(i) + '\t' + '\t'.join(str(x) for x in j) + '\n')

# to_write_profiles = []
# for i in unique_profiles_list:
#     i.insert(0, unique_profiles_list.index(i))
#     to_write_profiles.append('\t'.join(str(x) for x in i))
#
# with open(outputfile1, 'w') as outfile1:
#     outfile1.write('\n'.join(to_write_profiles))
#
df2 = pd.DataFrame.from_dict(combined_snp_positions, orient='columns')
df2.to_csv(path_or_buf="/home/conrad/Data/h_notebook_files/H6_treewas_maf005_genes_uncollapsed/gene_presence_absence_maf005_transposed_collapsed.txt", sep='\t')
