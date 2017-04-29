"""
for each high quality position:
    for each isolate:
        if an object with that position exists for the isolate:
            test if SNP, REF, or N
        elif an object with that position does not exist for the isolate:
            record REF

"""

from sys import argv

script, bwa6, novo6 ,last6 = argv

#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
# Define functions and classes
#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

# Define VcfVariants class. Allows easy access to different fields of a variant:
class VcfVariant(object):
    def __init__(self, chrom, pos, ref, alt, quality, dp, dp4, isolateid, variantid, caller):
        self.chrom = chrom
        self.pos = pos
        self.ref = ref
        self.alt = alt
        self.quality = quality
        self.dp = dp
        self.dp4 = dp4
        self.isolateid = isolateid
        self.variantid = variantid
        self.caller = caller
    def __repr__(self):
        return str(self.caller)

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
                temp1.append(vcf_file.split('/')[-1].split('_')[1].split('.')[0])
                variant_objects_list_1.append(VcfVariant(temp1[0], int(temp1[1]), temp1[2], temp1[3], float(temp1[4]), int(temp1[5]), temp1[6], temp1[7], temp1[8], temp1[9]))
    return variant_objects_list_1

# Define class for variant list objects:
class VariantListObject(object):
    def __init__(self, variant_list):
        self.variant_list = variant_list

    def positions_list(self):
        pos_list = []
        for i in self.variant_list:
            pos_list.append(i.pos)
        return pos_list
#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
bwa_objects = parse_vcf(bwa6)
novo_objects = parse_vcf(novo6)
last_objects = parse_vcf(last6)

bwa_list_object = VariantListObject(bwa_objects)
bwa_pos_list = bwa_list_object.positions_list()

novo_list_object = VariantListObject(novo_objects)
novo_pos_list = novo_list_object.positions_list()

last_list_object = VariantListObject(last_objects)
last_pos_list = last_list_object.positions_list()

bwa_intersect = [i for i in bwa_objects if i.pos in novo_pos_list and i.pos in last_pos_list]
print len(bwa_intersect)

novo_intersect = [i for i in novo_objects if i.pos in bwa_pos_list and i.pos in last_pos_list]
print len(novo_intersect)

last_intersect = [i for i in last_objects if i.pos in bwa_pos_list and i.pos in novo_pos_list]
print len(last_intersect)



#objects_intersect = [i for i in bwa_objects if i.pos in [j.pos for j in novo_objects] and i.pos in [k.pos for k in last_objects]]
#print len(objects_intersect)
