from sys import argv
import itertools

script, bwa6, novo6 ,last6 = argv

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

# Identify common positions between 3 aligners:
common_positions_25 = get_common_positions(bwa6, novo6, last6)

# Create lists of variant objects for each aligner:
bwa_objects = parse_vcf(bwa6)
novo_objects = parse_vcf(novo6)
last_objects = parse_vcf(last6)

# Identify variant objects that have positions in the common positions list created above:
bwa_filtered_list = get_common_variants(common_positions_25, bwa_objects)
novo_filtered_list = get_common_variants(common_positions_25, novo_objects)
last_filtered_list = get_common_variants(common_positions_25, last_objects)

# Call SNPs for each set of 3 objects in filtered lists:
dictt = {}
for bwasnp, novosnp, lastsnp in itertools.izip(bwa_filtered_list, novo_filtered_list, last_filtered_list):
    dictt[bwasnp.pos] = [bwasnp.ref, call_snps(bwasnp, novosnp, lastsnp)]

print dictt




#for i in bwa_filtered_list:
#    print str(i.pos) + '\t' + str(i.ref) + '\t' + str(i.alt)
