from sys import argv
from Bio import SeqIO

script, inputfasta = argv

dict1 = SeqIO.to_dict(SeqIO.parse(inputfasta, "fasta"))
print dict1
