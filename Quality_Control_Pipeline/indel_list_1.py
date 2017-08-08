from sys import argv
import itertools
import os.path

script, bwaindels, outputfile = argv

# Define VcfVariants class. Allows easy access to different fields of a variant:
class VcfVariant(object):
    def __init__(self, chrom, pos, ref, alt, quality, dp):
        self.chrom = chrom
        self.pos = pos
        self.ref = ref
        self.alt = alt
        self.quality = quality
        self.dp = dp


# Define function to parse VCF files and create VcfVariant class objects for each variant:
def parse_vcf(vcf_file):
    variant_objects_list_1 = []
    with open(vcf_file, 'r') as infile:
        for line in infile:
            if not line.startswith("#"):
                temp1 = line.split('\t')
                object1 = VcfVariant(temp1[0], int(temp1[1]), temp1[3], temp1[4], temp1[5], temp1[7].split(';')[3])
                variant_objects_list_1.append(object1)
    return variant_objects_list_1

bwa_indels = parse_vcf(bwaindels)

print "# indels: ", len(bwa_indels)

with open(outputfile, 'w') as outfile:
    outfile.write("Position" + '\t' + "Reference" + '\t' + "Alternate" + '\n')
    to_write = []
    for i in bwa_indels:
        to_write.append(str(i.pos) + '\t' + i.ref + '\t' + i.alt)
    outfile.write('\n'.join(to_write))
