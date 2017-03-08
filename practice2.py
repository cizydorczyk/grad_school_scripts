import itertools
from sys import argv
import os.path

script, inputfile = argv

"""Set High Quality SNP (HQS) Criteria:"""
HQS_DEPTH = 20
HQS_QUALITY = 30
# minimum number of reads in foward and reverse alternate reads:
HQS_ALT_THRESHOLD = 3
# Minimum distance (bp) to call a variant from ends/indels:
HQS_MIN_DIST = 150
# Length of reference:
HQS_REF_LENGTH = 6264404
# Min % difference between alt forward and alt reverse reads:
HQS_MIN_DP_PER = 0.20

"""Set Low Quality SNP (LQS) Criteria:"""
LQS_DEPTH = 15
LQS_QUALITY = 25
# minimum number of reads in forward and reverse alternate reads:
LQS_ALT_THRESHOLD = 3
# Minimum distance (bp) to call a variant from ends/indels:
LQS_MIN_DIST = 150
# Length of reference:
LQS_REF_LENGTH = 6264044
# Min % (decimal value) between alt f/w and ref f/w in dp4:
LQS_MIN_DP_PER = 0.20
# Min % (decimal value) reads supporting SNP:
LQS_REF_CUTOFF = 0.80
# MAX % (decimal value) reads supporting REF:
LQS_SNP_CUTOFF = 0.20

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

# Define function to create list of positions with indels:
def indel_position_list(variantlist):
    indel_position_list = []
    for variant in variantlist:
        if len(variant.ref) > 1 or len(variant.alt) > 1:
            if variant.pos not in indel_position_list:
                indel_position_list.append(variant.pos)
    return indel_position_list

# Define function to create list of high quality variant positions using output from above two functions:
def hq_positions_list(variantobjectslist, indelpositionlist):
    hq_pos_count = 0
    for variant in variantobjectslist:
        if variant.quality >= HQS_QUALITY:
            if variant.dp >= HQS_DEPTH and int(variant.dp4[2]) >= HQS_ALT_THRESHOLD and int(variant.dp4[3]) >= HQS_ALT_THRESHOLD and (float(variant.dp4[2]) / float(variant.dp4[3])) > HQS_MIN_DP_PER and (float(variant.dp4[3]) / float(variant.dp4[2])) > HQS_MIN_DP_PER and variant.pos >= HQS_MIN_DIST and variant.pos <= (HQS_REF_LENGTH - HQS_MIN_DIST) and variant.pos not in itertools.chain.from_iterable([range(indel - HQS_MIN_DIST, indel + HQS_MIN_DIST + 1) for indel in indelpositionlist]):
                hqpositions.append(variant.pos)
                hq_ref_dict[variant.pos] = variant.ref
                hq_pos_count += 1
    print "High quality positions identified in file: " + '\t' + str(hq_pos_count)


hqpositions = []
hq_ref_dict = {}
with open(inputfile, 'r') as infile:
    for line in infile:

        variant_objects = parse_vcf(line.strip())
        indel_list = indel_position_list(variant_objects)

        hq_positions_list(variant_objects, indel_list)
