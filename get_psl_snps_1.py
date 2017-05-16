from sys import argv
from Bio import SeqIO
import os.path

script, input_snp_file_list, reference_fasta, output_positions_file, outputfasta, outputdir = argv

# Define function to get common positions from a set of vcf files and output it as a list:
def get_common_positions(snp_file):

    positions1 = []

    with open(snp_file, 'r') as infile2:
        for line in infile2:
            if not line.startswith("Position"):
                line_pos = int(line.strip().split('\t')[0])
                if line_pos > psl_end:
                    continue
                elif psl_start <= line_pos:
                    positions1.append(line_pos)
                else:
                    continue

    return positions1

psl_start = 2453667
psl_end = 2472409

snp_file_list = []
with open(input_snp_file_list, 'r') as infile1:
    for line in infile1:
        snp_file_list.append(line.strip())

# Get all snp positions from psl operon:
positions = set()

for i in snp_file_list:
    isolate_positions = get_common_positions(i)
    for i in isolate_positions:
        positions.add(i)

# Parse sequence from reference fasta file:
for seq_record in SeqIO.parse(reference_fasta, "fasta"):
    ref_seq = seq_record.seq

# Create dictionary with reference bases for each hq position:
reference_dict = {}
for position in positions:
    if position not in reference_dict:
        reference_dict[position] = ref_seq[position - 1]

# Create positions list:
positions_list = list(positions)
positions_list.sort()

# Write hq SNP positions and reference bases to file:
with open(output_positions_file, 'w') as outfile1:
    outfile1.write("Position" + '\t' + "Ref" + '\n')
    for i in sorted(reference_dict):
        outfile1.write(str(i) + '\t' + str(reference_dict[i]) + '\n')

# Create each isolate's sequence based on snp positions identified above and write to file:
for i in snp_file_list:
    snp_dict = {}
    with open(i, 'r') as infile1:
        for line in infile1:
            if not line.startswith("Position"):
                position = int(line.strip().split('\t')[0])
                if position > psl_end:
                    continue
                elif psl_start <= position:
                    if position not in snp_dict:
                        snp_dict[position] = [line.strip().split('\t')[1], line.strip().split('\t')[2]]
                else:
                    continue

                # if position not in snp_dict:
                #     snp_dict[position] = [line.strip().split('\t')[1], line.strip().split('\t')[2]]

    sequence = ''
    for position1 in positions_list:
        if position1 in snp_dict:
            sequence += snp_dict[position1][1]
        else:
            sequence += reference_dict[position1]

    # Write sequences and snps called for each isolate to files:

    if os.path.isfile(outputfasta):
        with open(outputfasta, 'a') as outfile1:
            outfile1.write('>' + i.split('/')[-1].split('_')[0] + '\n')
            outfile1.write(sequence + '\n')
    elif not os.path.isfile(outputfasta):
        with open(outputfasta, 'w') as outfile1:
            outfile1.write('>' + i.split('/')[-1].split('_')[0] + '\n')
            outfile1.write(sequence + '\n')

    with open(outputdir + str(i.split('/')[-1].split('_')[0] + '_snps.txt'), 'w') as outfile2:
        outfile2.write('Position' + '\t' + 'Reference' + '\t' + 'Base_called' + '\n')
        for key in sorted(snp_dict):
            outfile2.write(str(key) + '\t' + str(snp_dict[key][0]) + '\t' + str(snp_dict[key][1] + '\n'))