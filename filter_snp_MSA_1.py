from sys import argv
import pandas

script, fasta_alignment, positions_file, output_alignment, output_positions = argv

input_dict = {}
with open(fasta_alignment, 'r') as infile:
    for line in infile:
        if line.startswith(">"):
            input_dict[line.strip()] = list(next(infile).strip())
positions_list = []
with open(positions_file, 'r') as infile2:
    for line in infile2:
        positions_list.append(line.strip())

reference = input_dict[">PAO1"]

df1 = pandas.DataFrame.from_dict(data=input_dict, orient='index')
# Sort dataframe in ascending order by row #:
df1.sort_index(axis=0, inplace=True)
# Remove reference sequence (PAO1) from dataframe:
df2 = df1.drop([">PAO1"])
# Count number of unique bases in each column:
uniques = df2.apply(lambda x: x.nunique())
# Create list of index positions of columns to drop:
to_drop = list(uniques[uniques==1].index)
# Remove positions from positions_list at indices identifed in to_drop:
for index in reversed(to_drop):
    positions_list.pop(index)
# Remove indices from to_drop from reference sequence:
for index in reversed(to_drop):
    reference.pop(index)

# Drop columns where the number of unique bases is 1:
df3 = df2.drop(uniques[uniques==1].index, axis=1)
# Transpose dataframe:
df4 = df3.transpose()
# Convert dataframe to dicionary, where isolate header is key, and list of non-identical positions is value:
output_dict = df4.to_dict(orient='list')

print "Positions removed: " + str(len(to_drop))
for key in output_dict:
    print len(output_dict[key])


with open(output_positions, 'w') as outfile1:
    for position in positions_list:
        outfile1.write(str(position) + '\n')

with open(output_alignment, 'w') as outfile2:
    for key, value in output_dict.iteritems():
        outfile2.write(key + '\n')
        outfile2.write(''.join(value) + '\n')
    outfile2.write('>PAO1' + '\n')
    outfile2.write(''.join(reference))
