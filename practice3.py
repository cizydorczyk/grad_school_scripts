from sys import argv

script, input_fasta, output_fasta = argv

fasta_record = []
with open(input_fasta, 'r') as infile:
    for line in infile:
        fasta_record.append(line)

with open(output_fasta, 'a+') as outfile:
    for line in fasta_record:
        outfile.write(line)
