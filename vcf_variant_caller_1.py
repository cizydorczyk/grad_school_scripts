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


#Define VcfVariants class;
class VcfVariants(object):
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


# Define Parser class, which takes a list of file names (full paths); VcfVariant class is called on each individual
# variant in each file
class Parser(object):
    def __init__(self, filename):
        self.filename = filename

    def vcf_parser_multiple(self):
        variant_dict = {}
        f = open(self.filename, 'r')
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
                temp_list1.append(self.filename)
                var_list1.append(temp_list1)
            # Call VcfVariant class on each parsed line (item in temp_list1)
        for j in var_list1:
            temp1 = VcfVariants(j[0], j[1], j[2], j[3], j[4], j[5], j[6], j[7], j[8])
            variant_list.append(temp1)

        variant_dict[self.filename] = variant_list
        f.close()
        return variant_dict

a = Parser(filename)
b = a.vcf_parser_multiple()