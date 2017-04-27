import sys
import os.path
import itertools
from Bio import SeqIO

script, bwa30list, novo30list, last30list, reference_fasta, output_positions_file = sys.argv

bwa30_list = []
with open(bwa30list, 'r') as infile1:
    for line in infile1:
        bwa30_list.append(line.strip())

novo30_list = []
with open(novo30list, 'r') as infile2:
    for line in infile2:
        novo30_list.append(line.strip())

last30_list = []
with open(last30list, 'r') as infile3:
    for line in infile3:
        last30_list.append(line.strip())

hq_pos_combined = set()
hq_pos_dict = {}

def parse_vcf(bwafile, novofile, lastfile):
    bwa_positions = []
    novo_positions = []
    last_positions = []

    with open(bwafile, 'r') as infile:
        for line in infile:
            if not line.startswith("#"):
                bwa_positions.append(int(line.strip().split('\t')[1]))

    with open(novofile, 'r') as infile:
        for line in infile:
            if not line.startswith("#"):
                novo_positions.append(int(line.strip().split('\t')[1]))

    with open(lastfile, 'r') as infile:
        for line in infile:
            if not line.startswith("#"):
                last_positions.append(int(line.strip().split('\t')[1]))

    hq_bwa_positions = set(bwa_positions)
    hq_novo_positions = set(novo_positions)
    hq_last_positions = set(last_positions)
    hq_positions = (hq_bwa_positions & hq_novo_positions & hq_last_positions)

    return hq_positions

# Identify hq positions from overlap of all 3 aligners at min qual 30:
for bwa, novo, last in itertools.izip(bwa30_list, novo30_list, last30_list):
    common_positions = parse_vcf(bwa, novo, last)
    for i in common_positions:
        if i not in hq_pos_combined:
            hq_pos_combined.add(i)

# Parse sequence from reference fasta file:
for seq_record in SeqIO.parse(reference_fasta, "fasta"):
    ref_seq = seq_record.seq

# Create dictionary with reference bases for each hq position:
for position in hq_pos_combined:
    if position not in hq_pos_dict:
        hq_pos_dict[position] = ref_seq[position-1]

# Create hq positions list from hq positions set and sort:
hq_pos_list = list(hq_pos_combined)
hq_pos_list.sort()

## Now have a sorted list of hq positions: hq_pos_list
## And a dictionary with reference bases for each hq position: hq_pos_dict

# Write hq SNP positions and reference bases to file:
with open(output_positions_file, 'w') as outfile1:
    for i in sorted(hq_pos_dict):
        outfile1.write(str(i) + '\t' + str(hq_pos_dict[i]) + '\n')

# Identify possible SNPs in each isolate with positions in hq_pos_list and
# Quality at least 25:

# bwa25_list = []
# with open(bwa25list, 'r') as infile1:
#     for line in infile1:
#         bwa25_list.append(line.strip())
#
# novo25_list = []
# with open(novo25list, 'r') as infile2:
#     for line in infile2:
#         novo25_list.append(line.strip())
#
# last25_list = []
# with open(last25list, 'r') as infile3:
#     for line in infile3:
#         last25_list.append(line.strip())
#
# # Identify hq positions from overlap of all 3 aligners at min qual 30:
# for bwa, novo, last in itertools.izip(bwa25_list, novo25_list, last25_list):
#     common_positions_25 = parse_vcf(bwa, novo, last)
#     # Parse entire VCF
#
#     # If
