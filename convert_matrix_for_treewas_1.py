from sys import argv
import numpy as np
import pandas as pd
import itertools

script, inputtxt, outputprofiles, outputguide, min_num_occurrences, max_num_occurrences = argv

df = pd.read_table(inputtxt)
df = df.transpose()
df.drop("Total_Count", axis=0, inplace=True)
df.columns = df.iloc[0]
df.drop("Gene", axis=0, inplace=True)

df.to_csv(path_or_buf="/home/conrad/Data/primary_project_3/gwas/G82/" + min_num_occurrences + '-' + max_num_occurrences + "_from_roary.txt", sep='\t')

df_dict = df.to_dict(orient="list")
df_list = list(df_dict)

guide_dict = {}
unique_profiles_list = []

combined_snp_positions = {}

isolate_names = df.index

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


with open(outputprofiles, 'w') as outfile1:
    for i, j in itertools.izip(range(0, len(unique_profiles_list)), unique_profiles_list):
        outfile1.write(str(i) + '\t' + '\t'.join(str(x) for x in j) + '\n')

with open(outputguide, 'w') as outfile2:
    for key, value in guide_dict.iteritems():
        outfile2.write(str(key) + '\t' + '\t'.join(value) + '\n')

df2 = pd.DataFrame.from_dict(combined_snp_positions, orient='columns')
df2.set_index(isolate_names, inplace=True)
print df2
df2.to_csv(path_or_buf="/home/conrad/Data/primary_project_3/gwas/G82/" + min_num_occurrences + '-' + max_num_occurrences + "_from_roary_condensed.txt", sep='\t')
