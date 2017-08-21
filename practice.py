from sys import argv

script, input_vcf = argv

################################################################################
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
################################################################################

vcf_objects = parse_vcf(input_vcf)
called_objects = []

for snp in vcf_objects:
    #print snp.pos, snp.ref, snp.alt
    if ((float(snp.dp4[2])+float(snp.dp4[3])) / (float(snp.dp4[0]) + float(snp.dp4[1]) + float(snp.dp4[2]) + float(snp.dp4[3]))) >= 0.80:
        # Call SNP so change nothing about the snp record
        #print "\t" + snp.ref, snp.alt
        called_objects.append(snp)
    elif 0.20 <= ((float(snp.dp4[2])+float(snp.dp4[3])) / (float(snp.dp4[0]) + float(snp.dp4[1]) + float(snp.dp4[2]) + float(snp.dp4[3]))) < 0.80:
        # Call N
        snp.alt = "N"
        #print "\t" + snp.ref, snp.alt
        called_objects.append(snp)
    elif ((float(snp.dp4[2])+float(snp.dp4[3])) / (float(snp.dp4[0]) + float(snp.dp4[1]) + float(snp.dp4[2]) + float(snp.dp4[3]))) < 0.20:
        # Call REF so change alt base to ref base
        snp.alt = snp.ref
        #print "\t" + snp.ref, snp.alt
        called_objects.append(snp)

for snp_ in called_objects:
    if snp_.alt == "N":
        print snp_.pos, snp_.ref, snp_.alt
