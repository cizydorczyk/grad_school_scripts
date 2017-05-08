from sys import argv

script, tab_annotation, GTF_intron_annotation, snps_file, ann_file_out = argv

codons = {'GCA':'A', 'GCC':'A', 'GCG':'A', 'GCT':'A', 'TGC':'C', 'TGT':'C', 'GAC':'D', 'GAT':'D',\
'GAA':'E', 'GAG':'E', 'TTC':'F', 'TTT':'F', 'GGA':'G', 'GGC':'G', 'GGG':'G', 'GGT':'G', 'CAC':'H', \
'CAT':'H', 'ATA':'I', 'ATC':'I', 'ATT':'I', 'AAA':'K', 'AAG':'K', 'TTA':'L', 'TTG':'L', 'CTA':'L', \
'CTC':'L', 'CTG':'L', 'CTT':'L', 'ATG':'M', 'AAC':'N', 'AAT':'N', 'CCA':'P', 'CCC':'P', 'CCG':'P', \
'CCT':'P', 'CAA':'Q', 'CAG':'Q', 'CGA':'R', 'CGC':'R', 'CGG':'R', 'CGT':'R', 'AGA':'R', 'AGG':'R', \
'TCA':'S', 'TCC':'S', 'TCG':'S', 'TCT':'S', 'AGC':'S', 'AGT':'S', 'ACA':'T', 'ACC':'T', 'ACG':'T', \
'ACT':'T', 'GTA':'V', 'GTC':'V', 'GTG':'V', 'GTT':'V', 'TGG':'W', 'TAC':'Y', 'TAT':'Y', 'TAA':'*', \
'TAG':'*', 'TGA':'*'}


# Intergenic annotation format:
# type (intergenic), start, end, length of intergenic region, left gene, right gene
GTF_intron_ann_dict = {}
with open(GTF_intron_annotation, 'r') as infile1:
    for line in infile1:
        if 'intergenic' in line:
            line = line.strip().split('\t')
            new_line = [line[2], line[3], line[4], str(int(line[4])-int(line[3]) + 1), line[8].split(';')[3].split(' ')[2].strip('"'), line[8].split(';')[4].split(' ')[2].strip('"')]
            key = (int(line[3]), int(line[4]))
            GTF_intron_ann_dict[key] = '\t'.join(new_line)
# Genic annotation format:
# Sequence, Locus Tag, Feature Type, Start, End, Strand, Name, Product Name, Accession,
# GI, Length (nucleotides), MW (predicted), Length (amino acids), Nucleotide Sequence,
# Amino Acid Sequence

class tab_annotation_object(object):
    def __init__(self, sequence, locus_tag, feature_type, start, end, strand, name, product_name, \
    accession, GI, length_nuc, MW, length_aa, nuc_seq, aa_seq, record):
        self.sequence = sequence
        self.locus_tag = locus_tag
        self.feature_type = feature_type
        self.start = start
        self.end = end
        self.strand = strand
        self.name = name
        self.product_name = product_name
        self.accession = accession
        self.GI = GI
        self.length_nuc = length_nuc
        self.MW = MW
        self.length_aa = length_aa
        self.nuc_seq = nuc_seq
        self.aa_seq =  aa_seq
        self.record = record


TAB_ann_dict = {}
with open(tab_annotation, 'r') as infile2:
    for line in infile2:
        if not line.startswith("#") and not line.startswith("Sequence"):
            line = line.strip().split('\t')
            key = (int(line[3]), int(line[4]))
            line.pop(15)
            line.pop(14)
            line.pop(12)
            line.pop(8)

            if len(line) == 15:
                TAB_ann_dict[key] = tab_annotation_object(line[0], line[1], line[2], \
                line[3], line[4], line[5], line[6],line[7], line[8], line[9], line[10], \
                line[11], line[12], line[13], line[14], line)
            elif len(line) == 14:
                line.append('-')
                TAB_ann_dict[key] = tab_annotation_object(line[0], line[1], line[2], \
                line[3], line[4], line[5], line[6],line[7], line[8], line[9], line[10], \
                line[11], line[12], line[13], line[14], line)



snps_dict = {}
with open(snps_file, 'r') as infile3:
    for line in infile3:
        if not line.startswith("Position"):
            snps_dict[int(line.split('\t')[0])] = [line.strip().split('\t')[1], \
            line.strip().split('\t')[2]]

def lookup_bounds(bounds, value):
    for min_, max_ in bounds:
        # Don't have to test value < min because this is being run on a sorted
        # list of tuples, so a value < min would have already been accepted in the
        # previous tuple, and if it's not in any tuples, the elif statement
        # covers that and returns something in the same step. It's broken down
        # into two statements to speed things up,
        # so that if the value is larger than the max of the tuple, it doesn't
        # even bother testing the minimum (making two comparisons):
        if value > max_:
            continue
        elif min_ <= value:
            return min_, max_
        else:
            return

sorted_TAB_ann_bounds = sorted(TAB_ann_dict)
sorted_GTF_intron_ann_bounds = sorted(GTF_intron_ann_dict)

to_file = []
for i in sorted(snps_dict):
    bound = lookup_bounds(sorted_TAB_ann_bounds, i)
    if bound is not None:
        alt_base = snps_dict[i][1]
        snp_pos = int(i)
        sample_sequence = list(TAB_ann_dict[bound].nuc_seq)
        # Zero-based index of the SNP in the gene:
        snp_index = (snp_pos - int(TAB_ann_dict[bound].start))
        # Ref codon:
        codon = ''

        if (snp_index + 1) % 3 == 0:
            codon = ''.join(sample_sequence[(snp_index-2):(snp_index + 1)])
        elif (snp_index + 1) % 3 != 0 and (snp_index + 2) % 3 == 0:
            codon = ''.join(sample_sequence[(snp_index-1):(snp_index+2)])
        elif (snp_index + 1) % 3 != 0 and (snp_index + 3) % 3 == 0:
            codon = ''.join(sample_sequence[snp_index:(snp_index + 3)])
        #print codon

        sample_sequence[snp_index] = alt_base

        # Alt codon:
        codon_alt = ''
        # Position in codon:
        position_in_codon = ''

        if (snp_index + 1) % 3 == 0:
            codon_alt = ''.join(sample_sequence[(snp_index-2):(snp_index + 1)])
            position_in_codon = 3
        elif (snp_index + 1) % 3 != 0 and (snp_index + 2) % 3 == 0:
            codon_alt = ''.join(sample_sequence[(snp_index-1):(snp_index+2)])
            position_in_codon = 2
        elif (snp_index + 1) % 3 != 0 and (snp_index + 3) % 3 == 0:
            codon_alt = ''.join(sample_sequence[snp_index:(snp_index + 3)])
            position_in_codon = 1
        #print codon_alt
        #print position_in_codon

        # Ref and alt aa:
        ref_aa = codons[codon]
        alt_aa = ''
        if 'N' not in codon_alt:
            alt_aa = codons[codon_alt]
        elif 'N' in codon_alt:
            alt_aa = '-'

        # Type of codon/aa change:
        codon_change = ''

        if ref_aa == alt_aa:
            codon_change = "synonymous"
        elif ref_aa != alt_aa:
            if alt_aa == '*':
                codon_change = "stop"
            elif alt_aa == '-':
                codon_change = "ambiguous"
            else:
                codon_change = "non-synonymous"
        #print codon_change

        output_line = [str(i), snps_dict[i][0], snps_dict[i][1], str((snp_index + 1)), codon, codon_alt, \
        str(position_in_codon), ref_aa, alt_aa, codon_change] + TAB_ann_dict[bound].record[0:13]

        output_line2 = '\t'.join(output_line)
        to_file.append(output_line2)

        # RRL = regulatory region left (gene)
        # RRR = regulatory region right (gene)
        # RRN = regulatory region none
    elif bound is None:
        type_ = ''
        bound2 = lookup_bounds(sorted_GTF_intron_ann_bounds, i)
        if bound2 is not None:
            min_, max_ = bound2
            min_bound = min_ + 150
            max_bound = max_ - 150
            if min_ <= i <= min_bound:
                type_ = 'RRL'
            elif max_bound <= i <= max_:
                type_ = 'RRR'
            else:
                type_ = 'RRN'

            intergenic_index = int(i) - (min_ + 1)
            to_file2 = str(i) + '\t' + snps_dict[i][0] + '\t' + snps_dict[i][1] + '\t' + str(intergenic_index) + '\t' + GTF_intron_ann_dict[bound2] + '\t' + type_

            to_file.append(to_file2)

with open(ann_file_out, 'w') as outfile37:
    outfile37.write('## Genic annotation format:' + '\n' + '## snp position, reference base, \
alternate base, position in gene (1-based), ref codon, alt codon, position in codon, \
ref aa, alt aa, type, Sequence, Locus Tag, Feature Type, Start, End, Strand, Name, Product Name, Accession, \
GI, Length (nuc), MW (predicted), Length (aa), Nucleotide Sequence, Amino Acid Sequence' + '\n' + '## Intergenic \
annotation format:' + '\n' + '## position, ref base, alt base, position in region (1-based), intergenic, start, end, length of intergenic region, left gene, \
 right gene, type' + '\n' + '## RRL = Regulatory Region Left, RRN = Regulatory Region None, RRR = Regulatory Region Right' + '\n')
    outfile37.write('\n'.join(to_file))