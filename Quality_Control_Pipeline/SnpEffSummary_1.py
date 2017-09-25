from sys import argv
import fnmatch

script, snpeff_ann_annotation, snpeff_eff_annotation, snpeff_ann_vcffile, snpeff_eff_annotation_onelinepervar, snpeff_ann_annotation_onelinepervar = argv

# Create table of codons:
codons = {'GCA':'A', 'GCC':'A', 'GCG':'A', 'GCT':'A', 'TGC':'C', 'TGT':'C', 'GAC':'D', 'GAT':'D',\
'GAA':'E', 'GAG':'E', 'TTC':'F', 'TTT':'F', 'GGA':'G', 'GGC':'G', 'GGG':'G', 'GGT':'G', 'CAC':'H', \
'CAT':'H', 'ATA':'I', 'ATC':'I', 'ATT':'I', 'AAA':'K', 'AAG':'K', 'TTA':'L', 'TTG':'L', 'CTA':'L', \
'CTC':'L', 'CTG':'L', 'CTT':'L', 'ATG':'M', 'AAC':'N', 'AAT':'N', 'CCA':'P', 'CCC':'P', 'CCG':'P', \
'CCT':'P', 'CAA':'Q', 'CAG':'Q', 'CGA':'R', 'CGC':'R', 'CGG':'R', 'CGT':'R', 'AGA':'R', 'AGG':'R', \
'TCA':'S', 'TCC':'S', 'TCG':'S', 'TCT':'S', 'AGC':'S', 'AGT':'S', 'ACA':'T', 'ACC':'T', 'ACG':'T', \
'ACT':'T', 'GTA':'V', 'GTC':'V', 'GTG':'V', 'GTT':'V', 'TGG':'W', 'TAC':'Y', 'TAT':'Y', 'TAA':'*', \
'TAG':'*', 'TGA':'*'}

class EffAnnotation(object):
    def __init__(self, chrom, pos, id_, ref, alt, effect, impact, funclass, codon, aa, aa_len, gene, biotype, coding, trid, rank):
        self.chrom = chrom
        self.pos = pos
        self.id = id_
        self.ref = ref
        self.alt = alt
        self.effect = effect
        self.impact = impact
        self.funclass = funclass
        self.codon = codon
        self.aa = aa
        self.aa_len = aa_len
        self.gene = gene
        self.biotype = biotype
        self.coding = coding
        self.trid = trid
        self.rank = rank

class AnnAnnotation(object):
    def __init__(self, chrom, pos, id_, ref, alt, allele, effect, impact, gene, geneid, feature, featureid, biotype, rank, hgvsc, hgvsp, cdna_pos,
    cdna_len, cds_pos, cds_len, aa_pos, aa_len, errors):
        self.chrom = chrom
        self.pos = pos
        self.id_ = id_
        self.ref = ref
        self.alt = alt
        self.allele = allele
        self.effect = effect
        self.impact = impact
        self.gene = gene
        self.geneid = geneid
        self.feature = feature
        self.featureid = featureid
        self.biotype = biotype
        self.rank = rank
        self.hgvsc = hgvsc
        self.hgvsp = hgvsp
        self.cdna_pos = cdna_pos
        self.cdna_len = cdna_len
        self.cds_pos = cds_pos
        self.cds_len = cds_len
        self.aa_pos = aa_pos
        self.aa_len = aa_len
        self.errors = errors

def parse_eff_annotation(eff_annotation_file):
    eff_objects = []
    with open(eff_annotation_file, 'r') as infile1:
        for line in infile1:
            if not line.startswith("CHROM\tPOS"):
                temp1 = line.strip().split('\t')
                eff_objects.append(EffAnnotation(temp1[0], temp1[1], temp1[2], temp1[3], temp1[4], temp1[5], temp1[6], temp1[7], temp1[8], temp1[9],
                temp1[10], temp1[11], temp1[12], temp1[13], temp1[14], temp1[15]))
    return eff_objects

def parse_ann_annotation(ann_annotation_file):
    ann_objects = []
    with open(ann_annotation_file, 'r') as infile2:
        for line in infile2:
            if not line.startswith("CHROM\tPOS"):
                temp1 = line.strip().split('\t')
                ann_objects.append(AnnAnnotation(temp1[0], temp1[1], temp1[2], temp1[3], temp1[4], temp1[5], temp1[6], temp1[7], temp1[8], temp1[9], temp1[10],
                temp1[11], temp1[12], temp1[13], temp1[14], temp1[15], temp1[16], temp1[17], temp1[18], temp1[19], temp1[20], temp1[21],
                temp1[22]))
    return ann_objects

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

                annotation1 = temp1[7].strip().split(';')
                annotation2 = fnmatch.filter(annotation1, 'ANN=*')
                if len(annotation2) == 0:
                    annotation2 = ['NNN|up/downstream eff. only|NA|NA|NA|NA|NA|NA|NA|NA|NA|NA|NA|NA|NA|NA']
                annotation3 = annotation2[0].split(',')
                annotation_dict = {}
                for i in annotation3:
                    temp_ann1 = i.split('|')[0:]
                    temp_ann2 = ANN_field(temp_ann1[0][-1], temp_ann1[1], temp_ann1[2], temp_ann1[3], temp_ann1[4], temp_ann1[5],
                    temp_ann1[6], temp_ann1[7], temp_ann1[8], temp_ann1[9], temp_ann1[10], temp_ann1[11], temp_ann1[12],
                    temp_ann1[13], temp_ann1[14], temp_ann1[15])
                    annotation_dict[i.split('|')[0][-1]] = temp_ann2

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

################################################################################
variant_objects = parse_vcf(snpeff_ann_vcffile)
header, isolate_genotype_order = parse_vcf_header(snpeff_ann_vcffile)

eff_objects = parse_eff_annotation(snpeff_eff_annotation)
ann_objects = parse_ann_annotation(snpeff_ann_annotation)
eff_oneline_objects = parse_eff_annotation(snpeff_eff_annotation_onelinepervar)
ann_oneline_objects = parse_ann_annotation(snpeff_ann_annotation_onelinepervar)

# # Number of SNPs:
number_of_variants = len(eff_objects)

# Number of variants with at least one 'N':
total_no_of_N_variants = 0
uniN = 0
biN = 0
triN = 0
tetraN= 0

test1 = []
for i in eff_objects:
    if 'N' in i.alt and len(i.alt.split(',')) == 1:
        uniN += 1
    elif 'N' in i.alt and len(i.alt.split(',')) == 2:
        biN += 1
    elif 'N' in i.alt and len(i.alt.split(',')) == 3:
        triN += 1
    elif 'N' in i.alt and len(i.alt.split(',')) == 4:
        tetraN += 1

total_no_of_N_variants = uniN + biN + triN + tetraN

# Number of bi and multiallelic SNPs:
biallelic_snps = 0
triallelic_snps = 0
tetraellelic_snps = 0

for i in eff_objects:
    if 'N' not in i.alt:
        if len(i.alt.split(',')) == 1:
            biallelic_snps += 1
        elif len(i.alt.split(',')) == 2:
            triallelic_snps += 1
        elif len(i.alt.split(',')) == 3:
            tetraellelic_snps += 1

multiallelic_snps = triallelic_snps + tetraellelic_snps

# Number of singletons:
singletons = 0

for i in variant_objects:
    if 'N' not in i.alt:
        if i.genotypes.count('.:.') == 154:
            singletons += 1

# # SnpEff Annotations Summary
# variant_type_dict = {'N_variants_(same_as_uniN)':0, 'Up/Downstream_eff_only':0}
#
# for i in variant_objects:
#     if len(i.alt) == 1 and 'N' in i.alt:
#         variant_type_dict['N_variants_(same_as_uniN)'] += 1
#     else:
#         # print i.record
#         # print i.annotation
#         for j in i.alt:
#             if j != 'N':
#                 try:
#                     # annotation37 = i.annotation[j].annotation
#                     if i.annotation[j].annotation not in variant_type_dict:
#                         variant_type_dict[i.annotation[j].annotation] = 1
#                     else:
#                         variant_type_dict[i.annotation[j].annotation] += 1
#                 except KeyError:
#                     variant_type_dict['Up/Downstream_eff_only'] += 1

mutated_genes_dict = {}
counted_positions = set()
for obj in ann_oneline_objects:
    if obj.pos not in counted_positions:
        counted_positions.add(obj.pos)
        if obj.geneid not in mutated_genes_dict:
            mutated_genes_dict[obj.geneid] = 1
        elif obj.geneid in mutated_genes_dict:
            mutated_genes_dict[obj.geneid] += 1
    elif obj.pos in counted_positions:
        pass

num_variants_dict = {1:0, 2:0, 3:0, 4:0, 5:0, 6:0, 7:0, 8:0, 9:0,
10:0, 20:0, 30:0, 40:0, 50:0, 60:0, '>60':0}

for i in mutated_genes_dict:
    if mutated_genes_dict[i] == 1:
        num_variants_dict[1] += 1
    elif mutated_genes_dict[i] == 2:
        num_variants_dict[2] += 1
    elif mutated_genes_dict[i] == 3:
        num_variants_dict[3] += 1
    elif mutated_genes_dict[i] == 4:
        num_variants_dict[4] += 1
    elif mutated_genes_dict[i] == 5:
        num_variants_dict[5] += 1
    elif mutated_genes_dict[i] == 6:
        num_variants_dict[6] += 1
    elif mutated_genes_dict[i] == 7:
        num_variants_dict[7] += 1
    elif mutated_genes_dict[i] == 8:
        num_variants_dict[8] += 1
    elif mutated_genes_dict[i] == 9:
        num_variants_dict[9] += 1
    elif mutated_genes_dict[i] == 10:
        num_variants_dict[10] += 1
    elif 10 < mutated_genes_dict[i] <= 20:
        num_variants_dict[20] += 1
    elif 20 < mutated_genes_dict[i] <= 30:
        num_variants_dict[30] += 1
    elif 30 < mutated_genes_dict[i] <= 40:
        num_variants_dict[40] += 1
    elif 40 < mutated_genes_dict[i] <= 50:
        num_variants_dict[50] += 1
    elif 50 < mutated_genes_dict[i] <= 60:
        num_variants_dict[60] += 1
    elif 60 < mutated_genes_dict[i]:
        num_variants_dict['>60'] += 1


# To test if the number of variant positions as summed from mutated genes dict
# is equal to total number of variant positions (counting each mutated positions
# once, even if it has multiple alleles):
# sum_ = 0
# for i in sorted(mutated_genes_dict):
#     sum_ += mutated_genes_dict[i]
#     if mutated_genes_dict[i] > 1:
#         print i, mutated_genes_dict[i]
# print sum_


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
print "Number of singeltons (no 'N' singletons):", singletons
print "\n"
# print "######### SnpEff Annotations Summary #########"
# for i in sorted(variant_type_dict):
#     print "Number of " + i, variant_type_dict[i]
# print "\n"
print "######### Numbers of genes/intergenic regions with X mutations #########"
for i in sorted(num_variants_dict):
    print str(i) + " mutations:", num_variants_dict[i]
