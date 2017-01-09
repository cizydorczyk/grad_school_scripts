import itertools
from sys import argv

script, filename = argv

"""Set High Quality SNP (HQS) Criteria:"""
HQS_DEPTH = 20
HQS_QUALITY = 30
# minimum number of reads in foward and reverse alternate reads:
HQS_ALT_THRESHOLD = 3
# Max % (decimal value) of Alt reads Ref reads can number to:
HQS_REF_MAX = 0.15
# Minimum distance (bp) to call a variant from ends/indels:
HQS_MIN_DIST = 150
# Length of reference:
HQS_REF_LENGTH = 6264044
# Min % (decimal value) between alt f/w and ref f/w in dp4:
HQS_MIN_DP_PER = 0.20

"""Set Low Quality SNP (LQS) Criteria:"""
LQS_DEPTH = 15
LQS_QUALITY = 25
# minimum number of reads in forward and reverse alternate reads:
LQS_ALT_THRESHOLD = 3
# Max % (decimal value) of Alt reads Ref reads can number to:
LQS_REF_MAX = 0.15
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


#Define VcfVariants class. Allows easy access to different fields of a variant.
class VcfVariant(object):
    def __init__(self, isolateid, chrom, pos, ref, alt, quality, dp, dp4, filename):
        self.isolateid = isolateid
        self.chrom = chrom
        self.pos = pos
        self.ref = ref
        self.alt = alt
        self.quality = quality
        self.dp = dp
        self. dp4 = dp4
        self.filename = filename

# Define vcf_parser function, which takes the input vcf file and outputs a list of VcfVariants
# class objects:
def vcf_parser(self):
    f = open(filename, 'r')
    var_list1 = []
    variant_list = []
    for line in f:
        if not line.startswith("#"):
            temp_list1 = []
            temp_list2 = line.split()
            temp_list1.append(temp_list2[0] + temp_list2[1])
            temp_list1.append(temp_list2[0])
            temp_list1.append(int(temp_list2[1]))
            temp_list1.append(temp_list2[3])
            temp_list1.append(temp_list2[4])
            temp_list1.append(float(temp_list2[5]))
            temp_list3 = temp_list2[7].split(";")
            for item in temp_list3:
                if item.startswith("DP="):
                    if len(item) == 5:
                        temp_list1.append(int(
                            temp_list3[temp_list3.index(item)][-2] + temp_list3[temp_list3.index(item)][-1]))
                    elif len(item) == 4:
                        temp_list1.append(int(temp_list3[temp_list3.index(item)][-1]))
                    elif len(item) == 6:
                        temp_list1.append(int(temp_list3[temp_list3.index(item)][-3] + temp_list3[temp_list3.index(item)][-2] + temp_list3[temp_list3.index(item)][-1]))
            for item in temp_list3:
                if item.startswith("DP4"):
                    temp_list4 = item[4:].split(",")
                    temp_list5 = []
                    for items in temp_list4:
                        temp_list5.append(int(items))
                    temp_list1.append(temp_list5)
            temp_list1.append(filename)
            var_list1.append(temp_list1)
    # Call VcfVariant class on each parsed line (item in temp_list1)
    for j in var_list1:
        temp1 = VcfVariant(j[0], j[1], j[2], j[3], j[4], j[5], j[6], j[7], j[8])
        variant_list.append(temp1)
    f.close()
    return variant_list

variantlist = vcf_parser(filename)

# Create list of positions with indels
def indel_pos_list(variantlist):
    indel_pos_list = []
    for variant in variantlist:
        if len(variant.ref) > 1 or len(variant.alt) > 1:
            if variant.pos not in indel_pos_list:
                indel_pos_list.append(variant.pos)
    return indel_pos_list

indel_list = indel_pos_list(variantlist)

# Used to sort high quality SNPs by position in hq_pos_list() function:
def byPos_key(variant):
    return variant.pos

# Creates list of high quality SNP positions by filtering all SNPs using HQS criteria
# at top of this file:
def hq_pos_list(variantlist, indel_list):

    hq_list = []
    hq_pos = []
    temp1 = []
    for variant in variantlist:
        if variant.dp >= HQS_DEPTH and variant.quality >= HQS_QUALITY and variant.dp4[2] >= HQS_ALT_THRESHOLD and variant.dp4[3] >= HQS_ALT_THRESHOLD and (variant.dp4[0] + variant.dp4[1]) < HQS_REF_MAX*(variant.dp4[2] + variant.dp4[3]) and variant.pos > HQS_MIN_DIST and variant.pos < (HQS_REF_LENGTH - HQS_MIN_DIST) and float(variant.dp4[2])/variant.dp4[3] > HQS_MIN_DP_PER and float(variant.dp4[3])/variant.dp4[2] > HQS_MIN_DP_PER and len(variant.ref) == 1 and len(variant.alt) == 1 and variant.pos not in itertools.chain.from_iterable([range(indel-HQS_MIN_DIST, indel + HQS_MIN_DIST + 1) for indel in indel_list]) and variant.pos not in hq_pos:
            hq_pos.append(variant.pos)
            temp1.append(variant)
    hq_var = sorted(temp1, key=byPos_key)
    hq_pos.sort()
    hq_list.append(hq_pos)
    hq_list.append(hq_var)
    return hq_list

hqposlist = hq_pos_list(variantlist, indel_list)

# Based on high quality SNP positions, relax criteria for SNPs and identify any other SNPs with lower quality at THOSE positions:
def lq_var_list(variantlist, indel_list, hqposlist):
    lqvarlist = []

    for variant in variantlist:
        if variant.pos in hqposlist[0] and variant.dp >= LQS_DEPTH and variant.quality >= LQS_QUALITY and variant.dp4[2] >= LQS_ALT_THRESHOLD and variant.dp4[3] >= LQS_ALT_THRESHOLD and (variant.dp4[0] + variant.dp4[1]) < LQS_REF_MAX *(variant.dp4[2] + variant.dp4[3]) and variant.pos > LQS_MIN_DIST and variant.pos < (LQS_REF_LENGTH - LQS_MIN_DIST) and float(variant.dp4[2]) / variant.dp4[3] > LQS_MIN_DP_PER and float(variant.dp4[3]) / variant.dp4[2] > LQS_MIN_DP_PER and len(variant.ref) == 1 and len(variant.alt) == 1 and variant.pos not in itertools.chain.from_iterable([range(indel - LQS_MIN_DIST, indel + LQS_MIN_DIST + 1) for indel in indel_list]) and variant.pos not in lqvarlist:
            lqvarlist.append(variant)
    return lqvarlist


lqvarlist = lq_var_list(variantlist, indel_list, hqposlist)
