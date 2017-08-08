from sys import argv

script, input_fasta, isolate_list, output_directory = argv

isolatelist = []
with open(isolate_list, 'r') as infile1:
    for line in infile1:
        isolatelist.append(line.strip())

fasta_dict = {}
with open(input_fasta, 'r') as infile2:
    for line in infile2:
        if '>' in line:
            fasta_dict[line[1:].strip()] = next(infile2)

for i in isolatelist:
    with open(output_directory + i + '_snps.fasta', 'a+') as outfile1:
        outfile1.write('>' + i + '\n' + fasta_dict[i] + '\n')
