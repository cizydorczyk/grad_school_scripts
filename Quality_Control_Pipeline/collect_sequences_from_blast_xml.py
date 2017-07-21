from sys import argv
from Bio.Blast import NCBIXML
import itertools
from decimal import *

script, blastxml, indiv_fasta, comb_fasta = argv
isolate_number = blastxml.split("/")[-1].split("_")[0]
gene_name = blastxml.split("/")[-1].split("_")[-1].split(".")[0]
print isolate_number
print gene_name
result_handle = open(blastxml, 'r')

# blast_record = NCBIXML.read(result_handle)
# for i in blast_record.alignments:
#     for hsp in i.hsps:
#         print len(hsp.sbjct)

blast_records = NCBIXML.parse(result_handle)
list_records = list(blast_records)

sequences_list = []
header_list = []
for blast_record in list_records:
    for record in blast_record.alignments:
        header = ">" + record.title.split(" ")[0]
        for hsp in record.hsps:
            sequences_list.append(hsp.sbjct)
            header += ' |' + "gene: " + gene_name + '|' + "% identity: " + str(round((Decimal(hsp.identities)/Decimal(hsp.align_length)) * Decimal(100), 2)) \
            + '|' + "matches: " + str(hsp.identities) + '|' + "alignment length: " + str(hsp.align_length) + '|' + "gaps: " + str(hsp.gaps) + '|'
            header_list.append(header)
            print header
# with open(outputfasta, 'w+') as outfile:
#     for sequence, number in itertools.izip(sequences_list, numbers):
#         outfile.write('>' + isolate_number + '_' + gene_name + '_' + str(number) + '\n' + sequence + '\n')
with open(indiv_fasta, 'w+') as outfile1:
    for header, sequence in itertools.izip(header_list, sequences_list):
        outfile1.write(header + '\n' + sequence + '\n')

with open(comb_fasta, 'a+') as outfile2:
    outfile2.write(header_list[0] + '\n' + sequences_list[0] + '\n')
