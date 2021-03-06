from sys import argv
import itertools
import pandas as pd
from collections import defaultdict


script, inputfasta, inputpositions = argv

positions_list = []
with open(inputpositions, 'r') as infile2:
    for line in infile2:
        positions_list.append(line.strip())

fasta_sequences = {}
with open(inputfasta, 'r') as infile1:
    for line in infile1:
        if '>' in line:
            fasta_sequences[line.strip()] = next(infile1).strip()

reference_sequence = fasta_sequences['>reference']
fasta_sequences.pop('>reference')

snp_profiles = {}
for key in fasta_sequences:
    snp_profiles[key] = []
    for i, j in itertools.izip(fasta_sequences[key], reference_sequence):
        if i == j:
            snp_profiles[key].append('1')
        else:
            snp_profiles[key].append('2')

df1 = pd.DataFrame.from_dict(snp_profiles, orient='index')
df1.columns = positions_list

snp_profiles_2 = df1.to_dict(orient='list')
for i in snp_profiles_2:
    snp_profiles_2[i] = ''.join(snp_profiles_2[i])

inverted = defaultdict(list)
for k, v in snp_profiles_2.items():
    inverted[v].append(k)


new_dict = {v[0]: k for k,v in inverted.items()}
keys = {v[0]: v for _, v in inverted.items()}

from pprint import pprint
pprint(inverted)
pprint(new_dict)
pprint(keys)                                                                                                                                                                                  
