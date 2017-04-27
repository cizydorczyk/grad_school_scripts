from sys import argv
import itertools
script, bwa6, novo6, last6 = argv

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

# Create list of common, quality 25 positions from 3 aligners
positions = parse_vcf(bwa6, novo6, last6)

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

# Define function to prase VCF files and create VcfVariant class objects for each variant:
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

# Create vcfvariant objects for each vcf file:
bwa_objects = vcf_objects(bwa6)
novo_objects = vcf_objects(novo6)
last_objects = vcf_objects(last6)

# Confirm that each vcf file for each aligner has a different amount of SNPs
# print len(bwa_objects)
# print len(novo_objects)
# print len(last_objects)


# For each list of vcf objects (one list per aligner), filter it and create a new list that contains only variants with positions found
# in the high quality positions identified above:
bwa_filtered_list = []
for i in bwa_objects:
    if i.pos in positions:
        bwa_filtered_list.append(i)
novo_filtered_list = []
for i in novo_objects:
    if i.pos in positions:
        novo_filtered_list.append(i)
last_filtered_list = []
for i in last_objects:
    if i.pos in positions:
        last_filtered_list.append(i)


## Confirm that each list, after filtering, contains the same positions in the same order
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

    # elif 0.20 < ((int(bwasnp.dp4[0]) + int(bwasnp.dp4[1])) / (int(bwasnp.dp4[0]) + int(bwasnp.dp4[1]) + int(bwasnp.dp4[2]) + int(bwasnp.dp4[3]))) < 0.80:
    #     print "AMBIGUOUS!"
    # else:
    #     print "Positions do not match among the 3 aligners!"

## Confirm that the length of each list is the same and equal to the length of common positions:

# print bwa_filtered_list[-5].pos
# print novo_filtered_list[-5].pos
# print last_filtered_list[-5].pos

# print len(bwa_filtered_list)
# print len(novo_filtered_list)
# print len(last_filtered_list)
# print len(positions)

## Code seems to work properly...
