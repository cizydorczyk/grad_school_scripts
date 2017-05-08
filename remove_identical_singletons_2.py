from sys import argv
import pandas

script, fasta, positions_file, output_fasta = argv

# Create dictionary of fasta sequences, with header as key and sequence as value:
infile_dict = {}
with open(fasta, 'r') as infile1:
    for line in infile1:
        if line.startswith(">"):
            infile_dict[line.strip()] = list(next(infile1).strip())

# Create list of high quality positions:
positions_list = []
# Reconstruct reference sequence from positions_file:
reference = []

with open(positions_file, 'r') as infile2:
    for line in infile2:
        if not line.startswith("Position"):
            positions_list.append(line.strip().split('\t')[0])
            reference.append(line.strip().split('\t')[1])

# Convert reference list into a series and label positions with high quality positions list:
ref_series = pandas.Series(reference, index=positions_list)

# Turn sequence dictionary into dataframe, with high quality positions as the column names:
df1 = pandas.DataFrame.from_dict(infile_dict, orient='index')
df1.columns = positions_list

# Convert dataframe into dictionary, with column name as key and list of values in each column as value:
columns = df1.to_dict(orient='list')

# Identify columns with singletons and identical columns, and add their keys to singletons list:
singletons = []

for key in columns:
    if len(key_set) == 2 and 'N' not in key_set:
        if 1 in [columns[key].count(key1) for key1 in key_set]:
            singletons.append(key)
    elif len(key_set) == 1:
        singletons.append(key)
    else:
        pass

# Drop columns in singeltons list from dataframe as well as the reference sequence:
df1.drop(singletons, axis=1, inplace=True)
ref_series.drop(singletons, axis=0, inplace=True)

# Convert dataframe and reference sequence (a series) to a dictionary and list, respectively (for easier writing to file):
df1_dict = df1.transpose().to_dict(orient='list')
ref_list = list(ref_series)

# Write reference sequence and isolate sequences to file:
with open(output_fasta, 'w') as outfile:
    outfile.write('>PAO1' + '\n')
    outfile.write(''.join(ref_list) + '\n')
    for key in df1_dict:
        outfile.write(key + '\n' + ''.join(df1_dict[key]) + '\n')

# Print some statistics to screen on what was called (probably want to write these down...):
print "Number of singletons and identical columns identified: " + str(len(singletons))
print "Original number of columns: " + str(len(columns))
print "New number of columns: " + str(len(columns) - len(singletons))
