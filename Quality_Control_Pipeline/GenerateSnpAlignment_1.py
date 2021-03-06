# Takes whole genome fasta and outputs snp alignment fasta...

from sys import argv
import pandas as pd

script, inputfile, output_fasta = argv

with open(inputfile, 'r') as infile:
    list1 = list(infile)

dict1 = {}
for num, i in enumerate(list1):
    if i.startswith(">"):
        # print i,
        tick = 1
        sequence = ''
        j = list1[num+tick]
        while not j.startswith(">"):
            sequence += j.strip()
            tick += 1
            try:
                j = list1[num+tick]
            except IndexError:
                break
        dict1[i.strip()[1:]] = list(sequence)

# for i in dict1:
#     if 'N' in dict1[i]:
#         print 'N present', i

df1 = pd.DataFrame.from_dict(dict1, orient='index')
print df1

# Convert dataframe into dictionary, with column name as key and list of values in each column as value:
columns = df1.to_dict(orient='list')

# Columns with N's (any #) and one other base:
n_tons = []

# Singletons:
singletons = []

# N's and at least two other bases:
n_var = []

# Identical columns:
identicals = []

# Biallelic SNPs
biallelic = 0

# Multiallelic SNPs
multiallelic =  []

for key in columns:
    key_set = set(columns[key])
    if len(key_set) == 1:
        identicals.append(key)
    elif len(key_set) == 2:
        if 'N' not in key_set:
            if 1 in [columns[key].count(key1) for key1 in key_set]:
                singletons.append(key)
            biallelic += 1
        elif 'N' in key_set:
            n_tons.append(key)
    elif len(key_set) > 2:
        if 'N' in key_set:
            n_var.append(key)
        else:
            # pass
            multiallelic.append(key)
        # print key_set

print 'original number of columns', len(columns)
print 'n_tons', len(n_tons)
print 'singletons', len(singletons)
print 'identicals', len(identicals)
print 'n_var', len(n_var)
print 'biallelic snps', str(biallelic)
print 'multiallelic snps', len(multiallelic)
print 'columns removed: identicals = ', len(identicals)

cols = list(df1)
nunique = df1.apply(pd.Series.nunique)
cols_to_drop = nunique[nunique == 1].index
df1.drop(cols_to_drop, axis=1, inplace=True)                                                                                                                                   

print df1

# Convert dataframe and reference sequence (a series) to a dictionary and list, respectively (for easier writing to file):
df1_dict = df1.transpose().to_dict(orient='list')
#ref_list = list(ref_series)

# Write reference sequence and isolate sequences to file:
with open(output_fasta, 'w') as outfile:
    #outfile.write('>PAO1' + '\n')
    #outfile.write(''.join(ref_list) + '\n')
    for key in df1_dict:
        outfile.write('>' + key + '\n' + ''.join(df1_dict[key]) + '\n')
