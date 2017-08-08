from sys import argv
import pandas

script, fasta, positions_file, output_fasta, output_positions, removed_positions_file = argv

# Create dictionary of fasta sequences, with header as key and sequence as value:
infile_dict = {}
with open(fasta, 'r') as infile1:
    for line in infile1:
        if line.startswith(">"):
            infile_dict[line.strip()] = list(next(infile1).strip())

# Create list of high quality positions:
positions_list = []

with open(positions_file, 'r') as infile2:
    for line in infile2:
        if not line.startswith("Position"):
            positions_list.append(line.strip().split('\t')[0])

#reference = infile_dict[">PAO1"]


# Convert reference list into a series and label positions with high quality positions list:
#ref_series = pandas.Series(reference, index=positions_list)

# Delete reference sequence from infile_dict:
#del infile_dict[">PAO1"]

# Turn sequence dictionary into dataframe, with high quality positions as the column names:
df1 = pandas.DataFrame.from_dict(infile_dict, orient='index')
df1.columns = positions_list

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
            pass
        multiallelic.append(key)
        print key_set

print 'original number of columns', len(columns)
print 'n_tons', len(n_tons)
print 'singletons', len(singletons)
print 'identicals', len(identicals)
print 'n_var', len(n_var)
print 'biallelic snps', str(biallelic)
print 'multiallelic snps', len(multiallelic)

# To remove any of the categories (n_tons, singletons, identicals, or n_var, just change what is being dropped:

df1.drop(n_tons, axis=1, inplace=True)
#ref_series.drop(n_tons, axis=0, inplace=True)

df1.drop(identicals, axis=1, inplace=True)
#ref_series.drop(identicals, axis=0, inplace=True)

df1.drop(singletons, axis=1, inplace=True)
df1.drop(multiallelic, axis=1, inplace=True)

for i in n_tons:
    positions_list.remove(i)

for i in identicals:
    positions_list.remove(i)

for i in singletons:
    positions_list.remove(i)

for i in multiallelic:
    positions_list.remove(i)

print 'total positions removed', len(n_tons) + len(identicals) + len(singletons) + len(multiallelic)
print 'new total number of columns', len(df1.columns)

# Convert dataframe and reference sequence (a series) to a dictionary and list, respectively (for easier writing to file):
df1_dict = df1.transpose().to_dict(orient='list')
#ref_list = list(ref_series)

# Write reference sequence and isolate sequences to file:
with open(output_fasta, 'w') as outfile:
    #outfile.write('>PAO1' + '\n')
    #outfile.write(''.join(ref_list) + '\n')
    for key in df1_dict:
        outfile.write(key + '\n' + ''.join(df1_dict[key]) + '\n')

# Write positions to new list:
with open(output_positions, 'w') as outfile2:
    for position in positions_list:
        outfile2.write(str(position) + '\n')


# Write removed positions to file:
with open(removed_positions_file, 'w') as outfile3:
    outfile3.write('\n'.join(n_tons) + '\n')
    outfile3.write('\n'.join(identicals))
