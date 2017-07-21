from sys import argv
from Bio import SeqIO

script, comb_fasta, output_fasta = argv
isolate_number = comb_fasta.split("/")[-1].split("_")[0]

fasta_records = SeqIO.parse(comb_fasta, "fasta")
fasta_list = list(fasta_records)

sequences_list = []
for sequence in fasta_list:
    sequences_list.append(sequence.seq)

sequence = ''
for i in sequences_list:
    sequence += i

with open(output_fasta, 'w+') as outfile:
    outfile.write(">" + isolate_number + "_" + "mlst_genes " + "|" + "acs" + "|"\
    + "aro" + "|" + "gua" + "|" + "mut" + "|" + "nuo" + "|" + "pps" + "|" + "trp"\
    + "\n" + str(sequence) + "\n")
