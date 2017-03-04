import itertools
from sys import argv
import os.path

script, log_file, vcf_file_list, positions_output_file, output_fasta = argv

"""Set High Quality SNP (HQS) Criteria:"""
HQS_DEPTH = 20
HQS_QUALITY = 30
# minimum number of reads in foward and reverse alternate reads:
HQS_ALT_THRESHOLD = 3
# Minimum distance (bp) to call a variant from ends/indels:
HQS_MIN_DIST = 150
# Length of reference:
HQS_REF_LENGTH = 6264404
# Min % difference between alt forward and alt reverse reads:
HQS_MIN_DP_PER = 0.20

"""Set Low Quality SNP (LQS) Criteria:"""
LQS_DEPTH = 15
LQS_QUALITY = 25
# minimum number of reads in forward and reverse alternate reads:
LQS_ALT_THRESHOLD = 3
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

# Define function to create list of positions with indels:
def indel_position_list(variantlist):
    indel_position_list = []
    for variant in variantlist:
        if len(variant.ref) > 1 or len(variant.alt) > 1:
            if variant.pos not in indel_position_list:
                indel_position_list.append(variant.pos)
    return indel_position_list

# Define function to create list of high quality variant positions using output from above two functions:
def hq_positions_list(variantobjectslist, indelpositionlist):
    hq_pos_count = 0
    for variant in variantobjectslist:
        if variant.quality >= HQS_QUALITY:
            if variant.dp >= HQS_DEPTH:
                if int(variant.dp4[2]) >= HQS_ALT_THRESHOLD and int(variant.dp4[3]) >= HQS_ALT_THRESHOLD:
                    if (float(variant.dp4[2]) / float(variant.dp4[3])) > HQS_MIN_DP_PER and (float(variant.dp4[3]) / float(variant.dp4[2])) > HQS_MIN_DP_PER:
                        if variant.pos >= HQS_MIN_DIST and variant.pos <= (HQS_REF_LENGTH - HQS_MIN_DIST):
                            if variant.pos not in itertools.chain.from_iterable([range(indel - HQS_MIN_DIST, indel + HQS_MIN_DIST + 1) for indel in indelpositionlist]):
                                if variant.pos not in set(hq_positions_list1):
                                    hq_positions_list1.append(variant.pos)
                                    hq_ref_dict[variant.pos] = variant.ref
                                    hq_pos_count += 1
    print "High quality positions identified in file: " + '\t' + str(hq_pos_count)
    infile0.write("High quality positions identified in file: " + '\t' + str(hq_pos_count) + '\n')

# Define function to create list of possible SNPs (relaxed filtering of variants; output is list of variant class objects):
def list_of_SNPs(variantobjectsset, indelpositionlist, hqpositionsset):
    list_of_snps = []
    for variant in variantobjectsset:
        if variant.pos in hqpositionsset:
            if variant.quality >= LQS_QUALITY:
                if variant.dp >= LQS_DEPTH:
                    if int(variant.dp4[2]) >= LQS_ALT_THRESHOLD and int(variant.dp4[3]) >= LQS_ALT_THRESHOLD:
                        if (float(variant.dp4[2]) / float(variant.dp4[3])) > LQS_MIN_DP_PER and (float(variant.dp4[3]) / float(variant.dp4[2])) > LQS_MIN_DP_PER:
                            if variant.pos >= LQS_MIN_DIST and variant.pos <= (LQS_REF_LENGTH - LQS_MIN_DIST):
                                if variant.pos not in itertools.chain.from_iterable([range(indel - LQS_MIN_DIST, indel + LQS_MIN_DIST + 1) for indel in indelpositionlist]):
                                    list_of_snps.append(variant)
    return list_of_snps

# Define function to call variants as SNPs or reference calls, or ambiguous:
def call_SNP(listofsnps, hqpositionslist, hqrefdict):
    sequence = ''
    snp_count = 0
    ref_count = 0
    amb_count = 0
    for position in hqpositionslist:
        for snp in listofsnps:
            if snp.pos == position:
                 if ((float(snp.dp4[2])+float(snp.dp4[3])) / (float(snp.dp4[0]) + float(snp.dp4[1]) + float(snp.dp4[2]) + float(snp.dp4[3]))) >= 0.80:
                     sequence += snp.alt
                     snp_count += 1
                     break
                 elif 0.20 <= ((float(snp.dp4[2])+float(snp.dp4[3])) / (float(snp.dp4[0]) + float(snp.dp4[1]) + float(snp.dp4[2]) + float(snp.dp4[3]))) < 0.80:
                     sequence += "N"
                     amb_count += 1
                     break
                 elif ((float(snp.dp4[2])+float(snp.dp4[3])) / (float(snp.dp4[0]) + float(snp.dp4[1]) + float(snp.dp4[2]) + float(snp.dp4[3]))) < 0.20:
                     sequence += snp.ref
                     ref_count += 1
                     break
        else:
            sequence += hqrefdict[position]
            ref_count += 1
    print "Alternate bases called: " + '\t' + str(snp_count)
    print "Reference bases called: " + '\t' + str(ref_count)
    print "Ambiguous bases called: " + '\t' + str(amb_count)

    infile0.write("Alternate bases called: " + '\t' + str(snp_count) + '\n')
    infile0.write("Reference bases called: " + '\t' + str(ref_count) + '\n')
    infile0.write("Ambiguous bases called: " + '\t' + str(amb_count) + '\n')


    return sequence

#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
# Create log file
infile0 = open(log_file, 'w')

# Initialize high quality positions list for all isolates:
hq_positions_list1 = []
# Initialize high quality positions and reference bases dictionary; used to call reference base at high
# quality positions that may be absent in a given vcf file:
hq_ref_dict = {}

# Parse all vcf files and identify high quality positions; append to hq_positions_list initialized above:
with open(vcf_file_list, 'r') as infile1:
    for line in infile1:
        print "Starting variant calling with file: " + '\t' + line.strip()
        infile0.write("Starting variant calling with file: " + line.strip() + '\n')
        print "Creating VcfVariant objects..."
        variant_objects = parse_vcf(line.strip())
        print "VcfVariant objects created: " + '\t' + str(len(variant_objects))
        infile0.write("VcfVariant objects created: " + '\t' + str(len(variant_objects)) + '\n')

        print "Creating list of indel positions..."
        indel_positions = indel_position_list(variant_objects)
        print "Indel positions identified: " + '\t' + str(len(indel_positions))
        infile0.write("Indel positions identified: " + '\t' + str(len(indel_positions)) + '\n')

        print "Creating list of high quality variant positions..."
        hq_positions_list(variant_objects, indel_positions) # doesn't need to return anything because this function appends to the global list "hq_positions_list1"...
        print "Total high quality variant positions identified: " + str(len(hq_positions_list1))
        infile0.write("Total high quality variant positions identified: " + '\t' + str(len(hq_positions_list1)) + '\n' + '\n')

# Sort hq_positions_list in ascending numerical order:
hq_positions_list1.sort()

# Write positions to output file:
with open(positions_output_file, 'w') as outfile0:
    for i in hq_positions_list1:
        outfile0.write(str(i) + '\n')
# Convert hq_positions_list to a set for faster membership testing in list_of_SNPs() function call:
hq_positions_set = set(hq_positions_list1)

# Parse each vcf file and identify possible SNPs using relaxed filtering criteria compared to hq_positions_list() function:
with open(vcf_file_list, 'r') as infile2:
    for line in infile2:
        print "Identifying possible SNPs in file: " + '\t' + line.strip()
        infile0.write("Identifying possible SNPs in file: " + '\t' + line.strip() + '\n')

        variant_objects = parse_vcf(line.strip())
        variant_objects_set = set(variant_objects)

        indel_positions = indel_position_list(variant_objects)

        print "Creating list of possible SNPs..."
        list_of_snps = list_of_SNPs(variant_objects_set, indel_positions, hq_positions_set)
        print "Possible SNPs identified: " + str(len(list_of_snps))
        infile0.write("Possible SNPs identified: " + '\t' + str(len(list_of_snps)) + '\n')

        print "Compiling isolate variant sequence..."
        sequence = call_SNP(list_of_snps, hq_positions_list1, hq_ref_dict)
        print "Bases in sequence: " + '\t' + str(len(sequence))
        infile0.write("Bases in sequence: " + '\t' + str(len(sequence)) + '\n' + '\n')

        if not os.path.isfile(output_fasta):
            with open(output_fasta, 'w') as outfile:
                outfile.write(">" + line.split('/')[-1].split('_')[0] + '\n')
                outfile.write(sequence + '\n')
        elif os.path.isfile(output_fasta):
            with open(output_fasta, 'a') as outfile:
                outfile.write(">" + line.split('/')[-1].split('_')[0] + '\n')
                outfile.write(sequence + '\n')



infile0.close()
#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
