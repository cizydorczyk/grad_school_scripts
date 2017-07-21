from sys import argv
import os.path
import itertools
from Bio import SeqIO

script, input_list_of_vcf_files, outputdirectory = argv

# Define VcfVariants class. Allows easy access to different fields of a variant:
class VcfVariant(object):
    def __init__(self, pos, ref, alt):
        self.pos = pos
        self.ref = ref
        self.alt = alt

# Define function to parse VCF files and create VcfVariant class objects for each variant:
def parse_vcf(vcf_file):
    variant_objects_list_1 = []
    with open(vcf_file, 'r') as infile:
        for line in infile:
            if not line.startswith("#"):
                temp1 = line.split('\t')
                variant_objects_list_1.append(VcfVariant(int(temp1[1]), temp1[3], temp1[4]))
    return variant_objects_list_1

vcf_file_list = []
with open(input_list_of_vcf_files, 'r') as infile1:
    for line in infile1:
        vcf_file_list.append(line.strip())

variant_dict = {}
for file_ in vcf_file_list:
    variant_list = parse_vcf(file_)
    isolate_snps = {}
    for i in variant_list:
        if i.pos not in variant_dict:
            variant_dict[i.pos] = i.ref
        isolate_snps[i.pos] = [i.ref, i.alt]
    with open(outputdirectory + file_.split("/")[-1].split("_")[0] + "_unannotated_snps.txt", 'w') as outfile3:
        outfile3.write("Position" + '\t' + "Reference" + '\t' + "Base_called" + '\n')
        for key in sorted(isolate_snps):
            outfile3.write(str(key) + '\t' + isolate_snps[key][0] + '\t' + isolate_snps[key][1] + '\n')
    #print file_
    #print len(variant_list)
print len(variant_dict)

with open(outputdirectory + 'all_snp_list.txt', 'w') as outfile1:
    for i in sorted(variant_dict):
        outfile1.write(str(i) + '\t' + variant_dict[i] + '\n')
