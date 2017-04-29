import sys
import os.path
import itertools
from Bio import SeqIO

script, bwa30list, novo30list, last30list, reference_fasta, output_positions_file, bwa25list, novo25list, last25list, outputdir, outputfasta = sys.argv

#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
# Define functions and classes
#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

# Define function to get common positions from a set of vcf files and output it as a list:
def get_common_positions(bwafile, novofile, lastfile):
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

# Define VcfVariants class. Allows easy access to different fields of a variant:
class VcfVariant(object):
    def __init__(self, chrom, pos, ref, alt, quality, dp, dp4, isolateid, variantid):
        self.chrom = chrom
        self.pos = pos
        self.ref = ref
        self.alt = alt
        self.quality = quality
        self.dp = dp
        self.dp4 = dp4
        self.isolateid = isolateid
        self.variantid = variantid

# Define function to parse VCF files and create VcfVariant class objects for each variant:
def parse_vcf(vcf_file):
    variant_objects_list_1 = []
    with open(vcf_file, 'r') as infile:
        for line in infile:
            if not line.startswith("#"):
                temp1 = line.split('\t')
                temp1.pop(9)
                temp1.pop(8)
                temp1.pop(6)
                temp1.pop(2)
                temp2 = temp1.pop(5).split(";")
                for item in temp2:
                    if item.startswith("DP="):
                        temp1.append(item[3:])
                    elif item.startswith("DP4="):
                        temp1.append(item[4:].split(','))
                temp1.append(vcf_file.split('/')[-1].split('_')[0])
                temp1.append(temp1[-1]+ "_" + temp1[1])
                variant_objects_list_1.append(VcfVariant(temp1[0], int(temp1[1]), temp1[2], temp1[3], float(temp1[4]), int(temp1[5]), temp1[6], temp1[7], temp1[8]))
    return variant_objects_list_1

# Define function to get common vcf variants, based on common positions list:
def get_common_variants(common_positions_list, variant_objects_list):
    common_list = []
    for i in variant_objects_list:
        if i.pos in common_positions_list:
            common_list.append(i)
    return common_list

# Define function to call snps:
def call_snps(bwasnp, novosnp, lastsnp):
    call = ''
    if ((float(bwasnp.dp4[0]) + float(bwasnp.dp4[1])) / (float(bwasnp.dp4[0]) + float(bwasnp.dp4[1]) + float(bwasnp.dp4[2]) + float(bwasnp.dp4[3]))) >= 0.80 \
    and ((float(novosnp.dp4[0]) + float(novosnp.dp4[1])) / (float(novosnp.dp4[0]) + float(novosnp.dp4[1]) + float(novosnp.dp4[2]) + float(novosnp.dp4[3]))) >= 0.80 \
    and ((float(lastsnp.dp4[0]) + float(lastsnp.dp4[1])) / (float(lastsnp.dp4[0]) + float(lastsnp.dp4[1]) + float(lastsnp.dp4[2]) + float(lastsnp.dp4[3]))) >= 0.80:
        call = bwasnp.ref
    elif ((float(bwasnp.dp4[0]) + float(bwasnp.dp4[1])) / (float(bwasnp.dp4[0]) + float(bwasnp.dp4[1]) + float(bwasnp.dp4[2]) + float(bwasnp.dp4[3]))) <= 0.20 \
    and ((float(novosnp.dp4[0]) + float(novosnp.dp4[1])) / (float(novosnp.dp4[0]) + float(novosnp.dp4[1]) + float(novosnp.dp4[2]) + float(novosnp.dp4[3]))) <= 0.20 \
    and ((float(lastsnp.dp4[0]) + float(lastsnp.dp4[1])) / (float(lastsnp.dp4[0]) + float(lastsnp.dp4[1]) + float(lastsnp.dp4[2]) + float(lastsnp.dp4[3]))) <= 0.20:
        call = bwasnp.alt
    else:
        call = 'N'
    return call

#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
# Create list of high SNP positions
#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

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

# Identify hq positions from overlap of all 3 aligners at min qual 30:
for bwa, novo, last in itertools.izip(bwa30_list, novo30_list, last30_list):
    common_positions = get_common_positions(bwa, novo, last)
    for i in common_positions:
        if i not in hq_pos_combined:
            hq_pos_combined.add(i)

# Parse sequence from reference fasta file:
for seq_record in SeqIO.parse(reference_fasta, "fasta"):
    ref_seq = seq_record.seq

# Create dictionary with reference bases for each hq position:
hq_pos_dict = {}
for position in hq_pos_combined:
    if position not in hq_pos_dict:
        hq_pos_dict[position] = ref_seq[position-1]

# Create hq positions list from hq positions set and sort:
hq_pos_list = list(hq_pos_combined)
hq_pos_list.sort()

# # Write hq SNP positions and reference bases to file:
# with open(output_positions_file, 'w') as outfile1:
#     outfile1.write("Position" + '\t' + "Ref" + '\n')
#     for i in sorted(hq_pos_dict):
#         outfile1.write(str(i) + '\t' + str(hq_pos_dict[i]) + '\n')

#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
# Identify SNPs in each isolate
#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

bwa25_list = []
with open(bwa25list, 'r') as infile1:
    for line in infile1:
        bwa25_list.append(line.strip())

novo25_list = []
with open(novo25list, 'r') as infile2:
    for line in infile2:
        novo25_list.append(line.strip())

last25_list = []
with open(last25list, 'r') as infile3:
    for line in infile3:
        last25_list.append(line.strip())

# For every isolate:
for bwa, novo, last in itertools.izip(bwa25_list, novo25_list, last25_list):

    # Identify common positions between 3 aligners:
    common_positions_25 = get_common_positions(bwa, novo, last)

    # Create lists of variant objects for each aligner:
    bwa_objects = parse_vcf(bwa)
    novo_objects = parse_vcf(novo)
    last_objects = parse_vcf(last)

    # Identify variant objects that have positions in the common positions list created above:
    bwa_filtered_list = get_common_variants(common_positions_25, bwa_objects)
    novo_filtered_list = get_common_variants(common_positions_25, novo_objects)
    last_filtered_list = get_common_variants(common_positions_25, last_objects)

    # Call SNPs for each set of 3 objects in filtered lists:
    snps = {}
    for bwasnp, novosnp, lastsnp in itertools.izip(bwa_filtered_list, novo_filtered_list, last_filtered_list):
        snps[bwasnp.pos] = [bwasnp.ref, call_snps(bwasnp, novosnp, lastsnp)]

    # Generate sequence based on high quality positions:
    sequence = ''
    for position in hq_pos_list:
        if position in snps:
            sequence += snps[position][1]
        else:
            sequence += hq_pos_dict[position]

    if os.path.isfile(outputfasta):
        with open(outputfasta, 'a') as outfile1:
            outfile1.write('>' + bwa.split('/')[-1].split('_')[0] + '\n')
            outfile1.write(sequence + '\n')
    elif not os.path.isfile(outputfasta):
        with open(outputfasta, 'w') as outfile1:
            outfile1.write('>' + bwa.split('/')[-1].split('_')[0] + '\n')
            outfile1.write(sequence + '\n')
    # Write sequences and snps called for each isolate to files:
    with open(outputdir+str(bwa.split('/')[-1].split('_')[0] +'_snps.txt'), 'w') as outfile2:
        outfile2.write('Position' + '\t' + 'Reference' + '\t' + 'Base_called' + '\n')
        for key in sorted(snps):
            outfile2.write(str(key) + '\t' + str(snps[key][0]) + '\t' + str(snps[key][1] + '\n'))







    # if os.path.isfile(outputfasta):
    #     with open(outputfasta, 'a') as outfile1:
    #         for key in sorted(snps):
    #             outfile1.write(str(key) + '\t' + str(snps[key][0]) + '\t' + str(snps[key][1] + '\n'))
    # elif not os.path.isfile(outputfasta):
    #     with open(outputfasta, 'w') as outfile1:
    #         outfile1.write('Position' + '\t' + 'Reference' + '\t' + 'Base_called' + '\n')
    #         for key in sorted(snps):
    #             outfile1.write(str(key) + '\t' + str(snps[key][0]) + '\t' + str(snps[key][1] + '\n'))

    # Write sequences and snps called for each isolate to files:
    # with open(outputdir+str(bwa.split('/')[-1].split('_')[0] +'snps.txt'), 'w') as outfile:
    #     outfile.write('Position' + '\t' + 'Reference' + '\t' + 'Base_called')
    #     to_write = []
    #     for key in sorted(snps):
    #         outfile.write(str(key) + '\t' + str(snps[key][0]) + '\t' + str(snps[key][1] + '\n')

    # if os.path.isfile(outputfasta):
    #     with open(outputfasta, 'a') as outfile1:
    #         outfile1.write('>' + bwa.split('/')[-1].split('_')[0] + '\n')
    #         outfile1.write(sequence + '\n')
    # elif not os.path.isfile(outputfasta):
    #     with open(outputfaasta, 'w') as outfile1:
    #         outfile1.write('>' + bwa.split('/')[-1].split('_')[0] + '\n')
    #         outfile1.write(sequence + '\n')
