import sys
import os.path
import itertools
from Bio import SeqIO

script, bwa30list, novo30list, last30list, reference_fasta, output_positions_file = sys.argv

#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
# Define functions and classes
#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

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
def vcf_objects(vcf_file):
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

# Define parse_vcf:
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
hq_pos_dict = {}

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
        outfile1.write("Position" + '\t' + "Ref" + '\n')
        outfile1.write(str(i) + '\t' + str(hq_pos_dict[i]) + '\n')
#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
# Identify SNPs in each isolate
#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

# Identify possible SNPs in each isolate with positions in hq_pos_list and
# Quality at least 25:

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
#
# # Identify hq positions from overlap of all 3 aligners at min qual 25:
for bwa, novo, last in itertools.izip(bwa25_list, novo25_list, last25_list):
    common_positions_25 = parse_vcf(bwa, novo, last)
    # Create vcfvariant objects for each vcf file:
    bwa_objects = vcf_objects(bwa)
    novo_objects = vcf_objects(novo)
    last_objects = vcf_objects(last)
    # For each list of vcf objects (one list per aligner), filter it and create a new list that contains only variants with positions found
    # in the high quality positions identified above:
    bwa_filtered_list = []
    for i in bwa_objects:
        if i.pos in common_positions_25:
            bwa_filtered_list.append(i)
    novo_filtered_list = []
    for i in novo_objects:
        if i.pos in common_positions_25:
            novo_filtered_list.append(i)
    last_filtered_list = []
    for i in last_objects:
        if i.pos in common_positions_25:
            last_filtered_list.append(i)

    for bwasnp, novosnp, lastsnp in itertools.izip(bwa_filtered_list, novo_filtered_list, last_filtered_list):
        #print (int(bwasnp.dp4[0]) + int(bwasnp.dp4[1]) + int(bwasnp.dp4[2]) + int(bwasnp.dp4[3]))

        if ((float(bwasnp.dp4[0]) + float(bwasnp.dp4[1])) / (float(bwasnp.dp4[0]) + float(bwasnp.dp4[1]) + float(bwasnp.dp4[2]) + float(bwasnp.dp4[3]))) >= 0.80 \
        and ((float(novosnp.dp4[0]) + float(novosnp.dp4[1])) / (float(novosnp.dp4[0]) + float(novosnp.dp4[1]) + float(novosnp.dp4[2]) + float(novosnp.dp4[3]))) >= 0.80 \
        and ((float(lastsnp.dp4[0]) + float(lastsnp.dp4[1])) / (float(lastsnp.dp4[0]) + float(lastsnp.dp4[1]) + float(lastsnp.dp4[2]) + float(lastsnp.dp4[3]))) >= 0.80:
            print "REF!" + '\t' + str(bwasnp.pos)
        elif ((float(bwasnp.dp4[0]) + float(bwasnp.dp4[1])) / (float(bwasnp.dp4[0]) + float(bwasnp.dp4[1]) + float(bwasnp.dp4[2]) + float(bwasnp.dp4[3]))) <= 0.20 \
        and ((float(novosnp.dp4[0]) + float(novosnp.dp4[1])) / (float(novosnp.dp4[0]) + float(novosnp.dp4[1]) + float(novosnp.dp4[2]) + float(novosnp.dp4[3]))) <= 0.20 \
        and ((float(lastsnp.dp4[0]) + float(lastsnp.dp4[1])) / (float(lastsnp.dp4[0]) + float(lastsnp.dp4[1]) + float(lastsnp.dp4[2]) + float(lastsnp.dp4[3]))) <= 0.20:
            print "SNP!"+ '\t' + str(bwasnp.pos)
        else:
            print "AMBIGUOUS!"+ '\t' + str(bwasnp.pos)
