# Takes whole genome fasta and outputs snp alignment fasta...

from sys import argv
import pandas as pd

script, inputfile, positions_file, outputtsv, outputfasta = argv

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

# Create list of high quality positions:
positions_list = []

with open(positions_file, 'r') as infile2:
    for line in infile2:
        if not line.startswith("Position"):
            positions_list.append(line.strip().split('\t')[0])

df1 = pd.DataFrame.from_dict(dict1, orient='index')

# Drop reference row from dataframe if exists:
try:
    df1.drop('reference', axis=0, inplace=True)
except KeyError:
    pass

# Define function to get number of each type of column:
def get_allelic(column_vector):
    key_set = set(column_vector)
    if len(key_set) is 1:
        if 'N' not in key_set:
            outcome = 'identical'
        elif 'N' in key_set:
            outcome = 'identical_N'
    elif len(key_set) is 2:
        if 'N' not in key_set:
            outcome = 'biallelic'
        elif 'N' in key_set:
            outcome = 'allele_+_N'
    elif len(key_set) is 3:
        if 'N' not in key_set:
            outcome = 'triallelic'
        elif 'N' in key_set:
            outcome = 'biallelic_+_N'
    elif len(key_set) is 4:
        if 'N' not in key_set:
            outcome = 'tetraallelic'
        elif 'N' in key_set:
            outcome = 'triallelic_+_N'
    elif len(key_set) is 5:
        if 'N' not in key_set:
            pass
        elif 'N' in key_set:
            outcome = 'tetraallelic_+_N'

    return outcome

# Apply function to get number of each type of column:
out = df1.apply(get_allelic, axis=0)
# Rename output from get_allelic function:
outnames = pd.Series(positions_list)
out.rename(index=outnames, inplace=True)

# Create list of column names to drop (change this to drop different things...):
to_drop_colnames = []
for i, j in out.iteritems():
    if j not in ['biallelic', 'triallelic', 'tetraallelic']:
        to_drop_colnames.append(i)

# Get number of each type of column:
allelic_dict = {}
for i, j in out.iteritems():
    if j not in allelic_dict:
        allelic_dict[j] = 1
    else:
        allelic_dict[j] += 1
for i in allelic_dict:
    print "Number of " + i + " sites:", allelic_dict[i]

# Define function to get number of singletons:
def get_singletons(column_vector):
    key_set = set(column_vector)
    if 'N' not in key_set:
        if len(key_set) is 2:
            if 1 in [list(column_vector).count(key1) for key1 in key_set]:
                outcome = 'singleton'
            else:
                outcome = 'not_singleton'
        else:
            outcome = 'not_singleton'
    elif 'N' in key_set:
        if len(key_set) is 2:
            if 1 in [list(column_vector).count(key1) for key1 in key_set]:
                outcome = 'n_singleton'
            else:
                outcome = 'n_not_singleton'
        else:
            outcome = 'n_not_singleton'

    return outcome

# Apply function to get number of singletons:
out2 = df1.apply(get_singletons, axis=0)
# Rename items in output from get_singletons function:
out2.rename(index=outnames, inplace=True)

# Count number of singletons
singletons_dict = {}
for i, j in out2.iteritems():
    if j not in singletons_dict:
        singletons_dict[j] = 1
    else:
        singletons_dict[j] += 1
for i in singletons_dict:
    print "Number of " + i + " sites:", singletons_dict[i]

# Name df columns and drop columns that are not biallelic SNPs:
df1.columns = positions_list
print "Number of original columns:", len(positions_list)
print "Number of columns to drop:", len(to_drop_colnames)
print "Number of columns remaining:", len(positions_list)-len(to_drop_colnames)
df1.drop(to_drop_colnames, axis=1, inplace=True)

# Transpose df1 to facilitate to_dict(); write dict as fasta:
df1_dict = df1.transpose().to_dict(orient='list')

with open(outputfasta, 'w') as outfile1:
    for key in df1_dict:
        outfile1.write(">" + key + '\n' + ''.join(df1_dict[key]) + '\n')

# Replace DNA bases with number values in the dataframe and create new df:
df1.replace(['A', 'C', 'G', 'T'], [1, 2, 3, 4], inplace=True)
print df1
# Write df to tsv:
df1.to_csv(path_or_buf=outputtsv, sep='\t')
