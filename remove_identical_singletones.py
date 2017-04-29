from sys import argv
import pandas

script, fasta, positions_file = argv

infile_dict = {}
with open(fasta, 'r') as infile1:
    for line in infile1:
        if line.startswith(">"):
            infile_dict[line.strip()] = list(next(infile1).strip())

positions_list = []
with open(positions_file, 'r') as infile2:
    for line in infile2:
        if not line.startswith("Position"):
            positions_list.append(line.strip().split('\t')[0])

df1 = pandas.DataFrame.from_dict(infile_dict, orient='index')
df1.columns = positions_list
