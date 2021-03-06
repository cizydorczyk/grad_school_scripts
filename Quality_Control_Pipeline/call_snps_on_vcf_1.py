from sys import argv

script, input_vcf, output_vcf = argv

################################################################################
# Define VcfVariants class. Allows easy access to different fields of a variant:
class VcfVariant(object):
    def __init__(self, chrom, pos, id_, ref, alt, qual, filter_, info, format_, extra, record, dp, dp4, mq):
        self.chrom = chrom
        self.pos = pos
        self.id_ = id_
        self.ref = ref
        self.alt = alt
        self.qual = qual
        self.filter = filter_
        self.info = info
        self.format_ = format_
        self.extra = extra
        self.record = record
        self.dp = dp
        self.dp4 = dp4
        self.mq = mq


# Define function to prase VCF files and create VcfVariant class objects for each variant:
def parse_vcf(vcf_file):
    variant_objects_list = []
    with open(vcf_file, 'r') as infile:
        for line in infile:
            if not line.startswith("#"):
                temp1 = line.strip().split('\t')
                temp3 = temp1[7].split(";")
                dp = ''
                dp4 = ''
                mq = ''
                for i in temp3:
                    if i.startswith("DP="):
                        dp = i[3:]
                    elif i.startswith("DP4="):
                        dp4 = i[4:].split(",")
                    elif i.startswith("MQ="):
                        mq = i[3:]

                variant_objects_list.append(VcfVariant(temp1[0], temp1[1], temp1[2],
                temp1[3], temp1[4], temp1[5], temp1[6], temp1[7], temp1[8], temp1[9], temp1,
                dp, dp4, mq))

    return variant_objects_list

def parse_vcf_header(vcf_input_file):
    header = ''
    with open(vcf_input_file, 'r') as infile2:
        for line in infile2:
            if line.startswith("#"):
                header += line

    return header

variant_objects = parse_vcf(input_vcf)
header = parse_vcf_header(input_vcf)

def call_snps(vcf_variant_objects_list):

    called_objects = []

    for snp in vcf_variant_objects_list:
        #print snp.pos, snp.ref, snp.alt
        if ((float(snp.dp4[2])+float(snp.dp4[3])) / (float(snp.dp4[0]) + float(snp.dp4[1]) + float(snp.dp4[2]) + float(snp.dp4[3]))) >= 0.90:
            # Call SNP so change nothing about the snp record
            called_objects.append(snp)
        elif 0.10 <= ((float(snp.dp4[2])+float(snp.dp4[3])) / (float(snp.dp4[0]) + float(snp.dp4[1]) + float(snp.dp4[2]) + float(snp.dp4[3]))) < 0.90:
            # Call N
            snp.alt = "N"
            snp.record[4] = "N"
            called_objects.append(snp)
        elif ((float(snp.dp4[2])+float(snp.dp4[3])) / (float(snp.dp4[0]) + float(snp.dp4[1]) + float(snp.dp4[2]) + float(snp.dp4[3]))) < 0.10:
            # Call REF so change alt base to ref base
            snp.alt = snp.ref
            snp.record[4] = snp.ref
            called_objects.append(snp)
    return called_objects

snps_called = call_snps(variant_objects)

for i in snps_called:
    if i.ref == i.alt:
        print i.ref, i.alt, i.dp4

with open(output_vcf, 'w') as outfile:
    outfile.write(header)
    to_write = []
    for snp in snps_called:
        to_write.append('\t'.join(snp.record))
    outfile.write('\n'.join(to_write))
