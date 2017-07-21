from sys import argv

script, input_fasta, output_fasta = argv
isolate_number = input_fasta.split("/")[-1].split("_")[0]

fasta_record = []
with open(input_fasta, 'r') as infile:
    for line in infile:
        fasta_record.append(line)

with open(output_fasta, 'a+') as outfile:
    outfile.write('>' + isolate_number + '\n')
    for line in fasta_record:
        if '>' not in line:
            outfile.write(line)
