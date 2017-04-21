from sys import argv

script, vcf250, vcf300 = argv

# Define VcfVariants class. Allows easy access to different fields of a variant:
class VcfVariant(object):
    def __init__(self, chrom, pos, ref, alt, info):
        self.chrom = chrom
        self.pos = pos
        self.ref = ref
        self.alt = alt
        self.info = info

# Define function to prase VCF files and create VcfVariant class objects for each variant:
def parse_vcf(vcf_file):
    variant_objects_list_1 = []
    with open(vcf_file, 'r') as infile:
        for line in infile:
            if not line.startswith('#'):
				temp1 = line.split('\t')
				temp1.pop(9)
				temp1.pop(8)
				temp1.pop(5)
				temp1.pop(2)
				temp1.pop(4)
				#print temp1
				variant_objects_list_1.append(VcfVariant(temp1[0], int(temp1[1]), temp1[2], temp1[3], temp1[4]))
	return variant_objects_list_1
				# temp1.pop(9)
                # temp1.pop(8)
                # temp1.pop(6)
                # temp1.pop(2)

                #temp2 = temp1.pop(5).split(";")
    #             for item in temp2:
    #                 if item.startswith("DP="):
    #                     temp1.append(item[3:])
    #                 elif item.startswith("DP4="):
    #                     temp1.append(item[4:].split(','))
    #             temp1.append(vcf_file.split('/')[-1].split('_')[0])
    #             temp1.append(temp1[-1]+ "_" + temp1[1])
    #             variant_objects_list_1.append(VcfVariant(temp1[0], int(temp1[1]), temp1[2], temp1[3], float(temp1[4]), int(temp1[5]), temp1[6], temp1[7], temp1[8]))
    # return variant_objects_list_1

print "parsing vcf 250"
vcf250 = parse_vcf(vcf250)
print len(vcf250)
print "parsing vcf 300"
vcf300 = parse_vcf(vcf300)
print len(vcf300)

two50_positions = []
three00_positions = []

for i in vcf250:
	two50_positions.append(i.pos)
for i in vcf300:
	three00_positions.append(i.pos)

set_250 = set(two50_positions)
set_300 = set(three00_positions)

for i in set_300:
	if i not in set_250:
		print "250 missing position " + str(i)

print "All 300 positions in 250 set"

for i in set_250:
	if i not in set_300:
		print "300 missing position " + str(i)
		
print "All 250 positions in 300 set"
