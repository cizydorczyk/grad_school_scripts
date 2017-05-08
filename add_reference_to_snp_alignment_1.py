from sys import argv

script, positions_file, raw_fasta = argv

reference = ''
with open(positions_file, 'r') as infile1:
    for line in infile1:
        if not line.startswith("Position"):
            reference += line.strip().split('\t')[1]

with open(raw_fasta, 'a') as outfile:
    outfile.write(">PAO1" + '\n' + reference)
