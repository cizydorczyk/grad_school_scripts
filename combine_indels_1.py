from sys import argv
import itertools
import os.path

script, bwaindels, novoindels, outputfile, logfile = argv

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
novo_indels = parse_vcf(novoindels)

print "# bwa indels: ", len(bwa_indels)
print "# novo indels: ", len(novo_indels)

# Define function to get common positions from a set of vcf files and output it as a list:
def get_common_positions(bwafile, novofile):
    bwa_positions = []
    novo_positions = []

    with open(bwafile, 'r') as infile:
        for line in infile:
            if not line.startswith("#"):
                bwa_positions.append(int(line.strip().split('\t')[1]))

    with open(novofile, 'r') as infile:
        for line in infile:
            if not line.startswith("#"):
                novo_positions.append(int(line.strip().split('\t')[1]))


    hq_bwa_positions = set(bwa_positions)
    hq_novo_positions = set(novo_positions)
    hq_positions = (hq_bwa_positions & hq_novo_positions)

    return hq_positions

common_indels = get_common_positions(bwaindels, novoindels)
print "# common indels: ", len(common_indels)

# Define function to get common vcf variants, based on common positions list:
def get_common_variants(common_positions_list, variant_objects_list):
    common_list = []
    for i in variant_objects_list:
        if i.pos in common_positions_list:
            common_list.append(i)
    return common_list

common_bwa_indels = get_common_variants(common_indels, bwa_indels)
common_novo_indels = get_common_variants(common_indels, novo_indels)

# Define function to identify identical indels, in terms of ref and alt (not qual/dp; position was done above):
def combine_indels(commonbwalist, commonnovolist):
    combined = []
    for bwa, novo in itertools.izip(commonbwalist, commonnovolist):
        if bwa.ref == novo.ref and bwa.alt == novo.alt:
            combined.append(bwa)
    return combined

combined_list = combine_indels(common_bwa_indels, common_novo_indels)
print "# identical indels: ", len(combined_list)

with open(outputfile, 'w') as outfile:
    outfile.write("Position" + '\t' + "Reference" + '\t' + "Alternate" + '\n')
    to_write = []
    for i in combined_list:
        to_write.append(str(i.pos) + '\t' + i.ref + '\t' + i.alt)
    outfile.write('\n'.join(to_write))

if os.path.isfile(logfile):
    with open(logfile, 'a') as outfile2:
        outfile2.write("# bwa indels:" + '\t' + str(len(bwa_indels)) + '\n' +
                       "# novo indels:" + '\t' + str(len(novo_indels)) + '\n' +
                       "# common indels based on position:" + '\t' + str(len(common_indels)) + '\n' +
                       "# identical indels based on ref and alt:" + '\t' + str(len(combined_list)) + '\n' + '\n')
elif not os.path.isfile(logfile):
    with open(logfile, 'w') as outfile2:
        outfile2.write("# bwa indels:" + '\t' + str(len(bwa_indels)) + '\n' +
                       "# novo indels:" + '\t' + str(len(novo_indels)) + '\n' +
                       "# common indels based on position:" + '\t' + str(len(common_indels)) + '\n' +
                       "# identical indels based on ref and alt:" + '\t' + str(len(combined_list)) + '\n' + '\n')