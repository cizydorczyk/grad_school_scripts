from sys import argv

script, input_positions, input_fasta, output_fasta_to_append = argv

with open(input_positions, 'r') as infile1:
    positions_list = list(infile1)

with open(input_fasta, 'r') as infile2:
    infile2.readline()
    pao1_seq_list = []
    for line in infile2:
        pao1_seq_list.append(line.strip())

pao1_seq = ''.join(pao1_seq_list)

aligned_seq = ''
for position in positions_list:
    aligned_seq += str(pao1_seq[int(position) -1])

with open(output_fasta_to_append, 'a') as outfile1:
    outfile1.write(">PAO1" + '\n')
    outfile1.write(aligned_seq)
