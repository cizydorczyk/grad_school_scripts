#### NOTE: There is nothing wrong with this script; it works perfectly. It's just difficult to further parse
# fields, such as getting the codons from the amino acid changes. That's why I switched to XXXXXXXXXXXX.
# This script works just fine for returning what it currently returns...

from sys import argv
import fnmatch

script, input_vcf, output_vcf = argv

################################################################################
# Define VcfVariants class. Allows easy access to different fields of a variant:
class VcfVariant(object):
    def __init__(self, chrom, pos, id_, ref, alt, qual, filter_, info, format_, extra, record, dp, dp4, mq, genotypes, annotation):
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
        self.genotypes = genotypes
        self.annotation = annotation

# class INFO_field(object):
#     def __init__(self, VDB, SGB, MQSB, MQ0F, MQ, RPB, MQB, BQB, DP, DP4, AN, AC, ANN, LOF):
#         self.VDB = VDB
#         self.SGB = SGB
#         self.MQSB = MQSB
#         self.MQ0F = MQ0F
#         self.MQ = MQ
#         self.RPB = RPB
#         self.MQB = MQB
#         self.BQB = BQB
#         self.DP = DP
#         self.DP4 = DP4
#         self.AN = AN
#         self.AC = AC
#         self.ANN = ANN
#         self.LOF = LOF

class ANN_field(object):
    def __init__(self, allele, annotation, putative_impact, gene_name, gene_id, feature_type, feature_id, transcript_biotype, rank_total, HGVSc, HGVSp,
    cDNApos_cDNAlen, CDSpos_CDSlen, PROTpos_PROTlen, dist_to_feature, errors):
        self.allele = allele
        self.annotation = annotation
        self.putative_impact = putative_impact
        self.gene_name = gene_name
        self.gene_id = gene_id
        self.feature_type = feature_type
        self.feature_id = feature_id
        self.transcript_biotype = transcript_biotype
        self.rank_total = rank_total
        self.HGVSc = HGVSc
        self.HGVSp = HGVSp
        self.cDNApos_cDNAlen = cDNApos_cDNAlen
        self.CDSpos_CDSlen = CDSpos_CDSlen
        self.PROTpos_PROTlen = PROTpos_PROTlen
        self.dist_to_feature = dist_to_feature
        self.errors = errors


# Define function to prase VCF files and create VcfVariant class objects for each variant:
def parse_vcf(vcf_file):
    variant_objects_list = []
    with open(vcf_file, 'r') as infile:
        for line in infile:
            if not line.startswith("#"):
                temp1 = line.strip().split('\t')
                alt = temp1[4].split(',')
                genotype_list = temp1[temp1.index('GT:PL') + 1:]

                # info_list = temp1[7].strip().split(';')
                # if len(info_list) == 13:
                #     info_list.append('NA')
                #     info = INFO_field(info_list[0], info_list[1], info_list[2], info_list[3], info_list[4], info_list[5],
                #     info_list[6], info_list[7], info_list[8], info_list[9], info_list[10], info_list[11], info_list[12],
                #     info_list[13])
                #
                #
                # elif len(info_list) == 14:
                #     info = INFO_field(info_list[0], info_list[1], info_list[2], info_list[3], info_list[4], info_list[5],
                #     info_list[6], info_list[7], info_list[8], info_list[9], info_list[10], info_list[11], info_list[12],
                #     info_list[13])

                annotation1 = temp1[7].strip().split(';')
                annotation2 = fnmatch.filter(annotation1, 'ANN=*')
                if len(annotation2) == 0:
                    annotation2 = ['NNN|up/downstream eff. only|NA|NA|NA|NA|NA|NA|NA|NA|NA|NA|NA|NA|NA|NA']
                annotation3 = annotation2[0].split(',')
                annotation_dict = {}
                for i in annotation3:
                    temp_ann1 = i.split('|')[0:]
                    # if len(temp_ann1) < 16:
                    #     print temp_ann1

                    temp_ann2 = ANN_field(temp_ann1[0][-1], temp_ann1[1], temp_ann1[2], temp_ann1[3], temp_ann1[4], temp_ann1[5],
                    temp_ann1[6], temp_ann1[7], temp_ann1[8], temp_ann1[9], temp_ann1[10], temp_ann1[11], temp_ann1[12],
                    temp_ann1[13], temp_ann1[14], temp_ann1[15])
                    annotation_dict[i.split('|')[0][-1]] = temp_ann2
                # print annotation_dict

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
                temp1[3], alt, temp1[5], temp1[6], temp1[7], temp1[8], temp1[9], temp1,
                dp, dp4, mq, genotype_list, annotation_dict))


    return variant_objects_list

def parse_vcf_header(vcf_input_file):
    header = ''
    isolate_genotype_order = []
    with open(vcf_input_file, 'r') as infile2:
        for line in infile2:
            if line.startswith("#"):
                header += line
                if line.startswith("#CHROM\tPOS"):
                    isolate_list = line.strip().split('\t')[9:]
                    for i in isolate_list:
                        isolate_genotype_order.append(i.strip().split('/')[-1].split('_')[0])

    return (header, isolate_genotype_order)

variant_objects = parse_vcf(input_vcf)
header, isolate_genotype_order = parse_vcf_header(input_vcf)

for i in variant_objects:
    if len(i.annotation) == 0:
        print i.record

# Number of SNPs:
number_of_variants = len(variant_objects)

# Number of variants with at least one 'N':
total_no_of_N_variants = 0
uniN = 0
biN = 0
triN = 0
tetraN= 0

for i in variant_objects:
    if 'N' in i.alt and len(i.alt) == 1:
        uniN += 1
    elif 'N' in i.alt and len(i.alt) == 2:
        biN += 1
    elif 'N' in i.alt and len(i.alt) == 3:
        triN += 1
    elif 'N' in i.alt and len(i.alt) == 4:
        tetraN += 1

total_no_of_N_variants = uniN + biN + triN + tetraN

# Number of bi and multiallelic SNPs:
biallelic_snps = 0
triallelic_snps = 0
tetraellelic_snps = 0

for i in variant_objects:
    if 'N' not in i.alt:
        if len(i.alt) == 1:
            biallelic_snps += 1
        elif len(i.alt) == 2:
            triallelic_snps += 1
        elif len(i.alt) == 3:
            tetraellelic_snps += 1

multiallelic_snps = triallelic_snps + tetraellelic_snps

# # Number of singletons:
# singletons = 0
#
# for i in variant_objects:
#     if 'N' not in i.alt:
#         if i.genotypes.count('.:.') == 154:
#             singletons += 1

# Number of non-/synonymous and intergenic SNPs:
variant_type_dict = {'N_variants_(same_as_uniN)':0, 'Up/Downstream_eff_only':0}

for i in variant_objects:
    if len(i.alt) == 1 and 'N' in i.alt:
        variant_type_dict['N_variants_(same_as_uniN)'] += 1
    else:
        # print i.record
        # print i.annotation
        for j in i.alt:
            if j != 'N':
                try:
                    # annotation37 = i.annotation[j].annotation
                    if i.annotation[j].annotation not in variant_type_dict:
                        variant_type_dict[i.annotation[j].annotation] = 1
                    else:
                        variant_type_dict[i.annotation[j].annotation] += 1
                except KeyError:
                    variant_type_dict['Up/Downstream_eff_only'] += 1

# Print summary to terminal:
print "Number of variants:", number_of_variants
print "\n"
print "Number of N variants (sum of uniN, biN, trin, & tetraN):", total_no_of_N_variants
print "Number of uniN variants (eg. Ref=A, Alt=N):", uniN
print "Number of biN variants (eg. Ref=A, Alt=T,N):", biN
print "Number of triN variants (eg. Ref=A, Alt=T,C,N):", triN
print "Number of tetraN variants (eg. Ref=A, Alt=T,C,G,N):", tetraN
print "\n"
print "Number of biallelic variants (eg. Ref=A, Alt=T):", biallelic_snps
print "Number of multiallelic variants (sum of tri & tetra allelics):", multiallelic_snps
print "Number of triallelic variants (eg. Ref=A, Alt=T,C):", triallelic_snps
print "Number of tetraellelic variants (eg. Ref=A, Alt=T,C,G):", tetraellelic_snps
print "\n"
print "Number of singeltons (npo 'N' singletons):", singletons
print "\n"
for i in sorted(variant_type_dict):
    print "Number of " + i, variant_type_dict[i]
