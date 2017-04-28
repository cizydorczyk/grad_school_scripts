from sys import argv
import numpy as np
import pandas as pd

script, bwafile, novofile, lastfile = argv

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
        return str(self.pos)

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

#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
bwa_objects = parse_vcf(bwafile)
novo_objects = parse_vcf(novofile)
last_objects = parse_vcf(lastfile)

variant_dict = {}

for i in bwa_objects:
    if i.pos not in variant_dict:
        if ((float(i.dp4[0]) + float(i.dp4[1])) / (float(i.dp4[0]) + float(i.dp4[1]) + float(i.dp4[2]) + float(i.dp4[3]))) >= 0.80:
            variant_dict[i.pos] = ['REF']
        elif ((float(i.dp4[0]) + float(i.dp4[1])) / (float(i.dp4[0]) + float(i.dp4[1]) + float(i.dp4[2]) + float(i.dp4[3]))) <= 0.20:
            variant_dict[i.pos] = ['SNP']
        else:
            variant_dict[i.pos] = ['N']
for i in novo_objects:
    if i.pos not in variant_dict:
        if ((float(i.dp4[0]) + float(i.dp4[1])) / (float(i.dp4[0]) + float(i.dp4[1]) + float(i.dp4[2]) + float(i.dp4[3]))) >= 0.80:
            variant_dict[i.pos] = [np.nan, 'REF']
        elif ((float(i.dp4[0]) + float(i.dp4[1])) / (float(i.dp4[0]) + float(i.dp4[1]) + float(i.dp4[2]) + float(i.dp4[3]))) <= 0.20:
            variant_dict[i.pos] = [np.nan, 'SNP']
        else:
            variant_dict[i.pos] = [np.nan, 'N']
    elif i.pos in variant_dict:
        if ((float(i.dp4[0]) + float(i.dp4[1])) / (float(i.dp4[0]) + float(i.dp4[1]) + float(i.dp4[2]) + float(i.dp4[3]))) >= 0.80:
            variant_dict[i.pos].append('REF')
        elif ((float(i.dp4[0]) + float(i.dp4[1])) / (float(i.dp4[0]) + float(i.dp4[1]) + float(i.dp4[2]) + float(i.dp4[3]))) <= 0.20:
            variant_dict[i.pos].append('SNP')
        else:
            variant_dict[i.pos].append('N')
for i in last_objects:
    if i.pos not in variant_dict:
        if ((float(i.dp4[0]) + float(i.dp4[1])) / (float(i.dp4[0]) + float(i.dp4[1]) + float(i.dp4[2]) + float(i.dp4[3]))) >= 0.80:
            variant_dict[i.pos] = [np.nan, 'REF']
        elif ((float(i.dp4[0]) + float(i.dp4[1])) / (float(i.dp4[0]) + float(i.dp4[1]) + float(i.dp4[2]) + float(i.dp4[3]))) <= 0.20:
            variant_dict[i.pos] = [np.nan, 'SNP']
        else:
            variant_dict[i.pos] = [np.nan, 'N']
    elif i.pos in variant_dict:
        if ((float(i.dp4[0]) + float(i.dp4[1])) / (float(i.dp4[0]) + float(i.dp4[1]) + float(i.dp4[2]) + float(i.dp4[3]))) >= 0.80:
            variant_dict[i.pos].append('REF')
        elif ((float(i.dp4[0]) + float(i.dp4[1])) / (float(i.dp4[0]) + float(i.dp4[1]) + float(i.dp4[2]) + float(i.dp4[3]))) <= 0.20:
            variant_dict[i.pos].append('SNP')
        else:
            variant_dict[i.pos].append('N')

for key in variant_dict:
    if len(variant_dict[key]) == 2:
        variant_dict[key].append(np.nan)
    elif len(variant_dict[key]) == 1:
        variant_dict[key].append(np.nan)
        variant_dict[key].append(np.nan)

vardf = pd.DataFrame.from_dict(variant_dict, orient='index')
vardf.sort_index(inplace=True)
vardf.columns = ['bwa', 'novo', 'last']
vardf.dropna(axis=0, inplace=True)


print vardf.iloc[i]
