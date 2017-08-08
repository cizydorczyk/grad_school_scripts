from sys import argv

script, tab_annotation, GTF_intron_annotation, indeltabfile, outfile = argv

# Intergenic annotation format:
# type (intergenic), start, end, length of intergenic region, left gene, right gene
# Parse intergenic annotation file:

class IntergenicObject(object):
    def __init__(self, type, start, end, length, left_gene, right_gene, record):
        self.type = type
        self.start = start
        self.end = end
        self.length = length
        self.left_gene = left_gene
        self.right_gene = right_gene
        self.record = record

GTF_intron_ann_dict = {}
with open(GTF_intron_annotation, 'r') as infile1:
    for line in infile1:
        if 'intergenic' in line:
            line = line.strip().split('\t')
            new_line = [line[2], line[3], line[4], str(int(line[4])-int(line[3]) + 1), line[8].split(';')[3].split(' ')[2].strip('"'),
                        line[8].split(';')[4].split(' ')[2].strip('"')]
            key = (int(line[3]), int(line[4]))
            GTF_intron_ann_dict[key] = IntergenicObject(line[2], line[3], line[4], str(int(line[4])-int(line[3]) + 1),
                                                        line[8].split(';')[3].split(' ')[2].strip('"'), line[8].split(';')[4].split(' ')[2].strip('"'),
                                                        '\t'.join(new_line))

# Genic annotation format:
# Sequence, Locus Tag, Feature Type, Start, End, Strand, Name, Product Name, Accession,
# GI, Length (nucleotides), MW (predicted), Length (amino acids), Nucleotide Sequence,
# Amino Acid Sequence

# Create class for gene annotation objects:
class tab_annotation_object(object):
    def __init__(self, sequence, locus_tag, feature_type, start, end, strand, name, product_name,
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
        self.aa_seq = aa_seq
        self.record = record

# Parse gene annotation file, creating class object for each annotation:
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

            TAB_ann_dict[key] = tab_annotation_object(line[0], line[1], line[2],
                line[3], line[4], line[5], line[6],line[7], line[8], line[9], line[10],
                line[11], line[12], line[13], line[14], line)

# Create dictionary of indels:
indels_dict = {}
with open(indeltabfile, 'r') as infile3:
    for line in infile3:
        if not line.startswith("Position"):
            position = int(line.split('\t')[0])
            if position not in indels_dict:
                indels_dict[position] = [line.split('\t')[1], line.strip().split('\t')[2]]

# Insertion or deletion:
annotation_dict = {}
for i in indels_dict:
    annotation_dict[i] = [str(i)]
    annotation_dict[i].append(indels_dict[i][0])
    annotation_dict[i].append(indels_dict[i][1])
    if len(indels_dict[i][0]) > len(indels_dict[i][1]):
        annotation_dict[i].append('deletion')
    else:
        annotation_dict[i].append('insertion')

# Length indel:
for i in indels_dict:
    if len(indels_dict[i][0]) > len(indels_dict[i][1]):
        annotation_dict[i].append(str(len(indels_dict[i][0])-1))
    else:
        annotation_dict[i].append(str(len(indels_dict[i][1])-1))

# Where indel occurs:

# Create function for testing if snp is in range of a gene/intergenic region:
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
# Create function to identify where an indel starts/stops, and specify whether it is intergenic or not:
def indel_location(indels_dict):

    sorted_TAB_ann_bounds = sorted(TAB_ann_dict)
    sorted_GTF_intron_ann_bounds = sorted(GTF_intron_ann_dict)

    for i in indels_dict:
        if 'insertion' in annotation_dict[i]: # Insertion
            alt_seq = indels_dict[i][1]
            alt_length = len(alt_seq) - 1

            alt_start = i + 1
            alt_stop = alt_start + alt_length - 1

            alt_start_bound = lookup_bounds(sorted_TAB_ann_bounds, alt_start)
            alt_stop_bound = lookup_bounds(sorted_TAB_ann_bounds, alt_stop)

            if alt_start_bound is not None:
                alt_start_region = TAB_ann_dict[alt_start_bound].locus_tag
                locus = 'genic'

                intergenic_type1 = 'CDS'

                min_, max_ = alt_start_bound
                position_in_gene = alt_start - min_ + 1

            else:
                alt_start_intergenic_bound = lookup_bounds(sorted_GTF_intron_ann_bounds, alt_start)
                alt_start_region = GTF_intron_ann_dict[alt_start_intergenic_bound].left_gene + ':' + GTF_intron_ann_dict[
                    alt_start_intergenic_bound].right_gene
                locus = 'intergenic'

                min_, max_ = alt_start_intergenic_bound
                min_bound = min_ + 150
                max_bound = max_ - 150
                if min_ <= i + 1 <= min_bound:
                    intergenic_type1 = 'RRL'
                elif max_bound <= i + 1 <= max_:
                    intergenic_type1 = 'RRR'
                else:
                    intergenic_type1 = 'RRN'

                position_in_gene = alt_start - min_ + 1

            if alt_stop_bound is not None:
                alt_stop_region = TAB_ann_dict[alt_stop_bound].locus_tag
                locus2 = 'genic'

                intergenic_type2 = 'CDS'

                min_, max_ = alt_stop_bound
                position_in_gene2 = alt_stop - min_ + 1

            else:
                alt_stop_intergenic_bound = lookup_bounds(sorted_GTF_intron_ann_bounds, alt_stop)
                alt_stop_region = GTF_intron_ann_dict[alt_stop_intergenic_bound].left_gene + ':' + GTF_intron_ann_dict[
                    alt_stop_intergenic_bound].right_gene
                locus2 = 'intergenic'

                min_, max_ = alt_stop_intergenic_bound
                min_bound = min_ + 150
                max_bound = max_ - 150
                if min_ <= i + 1 <= min_bound:
                    intergenic_type2 = 'RRL'
                elif max_bound <= i + 1 <= max_:
                    intergenic_type2 = 'RRR'
                else:
                    intergenic_type2 = 'RRN'

                position_in_gene2 = alt_stop - min_ + 1

            alt_ann = alt_start_region + '/' + alt_stop_region
            locus_annotation = locus + '/' + locus2
            type_ = intergenic_type1 + '/' + intergenic_type2
            position_ = str(position_in_gene) + '/' + str(position_in_gene2)

            annotation_dict[i].append(locus_annotation)
            annotation_dict[i].append(alt_ann)
            annotation_dict[i].append(type_)
            annotation_dict[i].append(position_)

        else: # Deletion
            ref_seq = indels_dict[i][0]
            ref_length = len(ref_seq) - 1

            ref_start = i + 1
            ref_stop = ref_start + ref_length - 1

            ref_start_bound = lookup_bounds(sorted_TAB_ann_bounds, ref_start)
            ref_stop_bound = lookup_bounds(sorted_TAB_ann_bounds, ref_stop)

            if ref_start_bound is not None:
                ref_start_region = TAB_ann_dict[ref_start_bound].locus_tag
                locus = 'genic'

                intergenic_type1 = 'CDS'

                min_, max_ = ref_start_bound
                position_in_gene = ref_start - min_ + 1

            else:
                ref_start_intergenic_bound = lookup_bounds(sorted_GTF_intron_ann_bounds, ref_start)
                ref_start_region = GTF_intron_ann_dict[ref_start_intergenic_bound].left_gene + ':' + GTF_intron_ann_dict[
                    ref_start_intergenic_bound].right_gene
                locus = 'intergenic'

                min_, max_ = ref_start_intergenic_bound
                min_bound = min_ + 150
                max_bound = max_ - 150
                if min_ <= i + 1 <= min_bound:
                    intergenic_type1 = 'RRL'
                elif max_bound <= i + 1 <= max_:
                    intergenic_type1 = 'RRR'
                else:
                    intergenic_type1 = 'RRN'

                position_in_gene = ref_start - min_ + 1

            if ref_stop_bound is not None:
                ref_stop_region = TAB_ann_dict[ref_stop_bound].locus_tag
                locus2 = 'genic'

                intergenic_type2 = 'CDS'

                min_, max_ = ref_stop_bound
                position_in_gene2 = ref_stop - min_ + 1

            else:
                ref_stop_intergenic_bound = lookup_bounds(sorted_GTF_intron_ann_bounds, ref_stop)
                ref_stop_region = GTF_intron_ann_dict[ref_stop_intergenic_bound].left_gene + ':' + GTF_intron_ann_dict[
                    ref_stop_intergenic_bound].right_gene
                locus2 = 'intergenic'

                min_, max_ = ref_stop_intergenic_bound
                min_bound = min_ + 150
                max_bound = max_ - 150
                if min_ <= i + 1 <= min_bound:
                    intergenic_type2 = 'RRL'
                elif max_bound <= i + 1 <= max_:
                    intergenic_type2 = 'RRR'
                else:
                    intergenic_type2 = 'RRN'

                position_in_gene2 = ref_stop - min_ + 1

            ref_ann = ref_start_region + '/' + ref_stop_region
            locus_annotation = locus + '/' + locus2
            type_ = intergenic_type1 + '/' + intergenic_type2
            position_ = str(position_in_gene) + '/' + str(position_in_gene2)

            annotation_dict[i].append(locus_annotation)
            annotation_dict[i].append(ref_ann)
            annotation_dict[i].append(type_)
            annotation_dict[i].append(position_)

indel_location(indels_dict)

# Frameshift or inframe:

for i in annotation_dict:
    if int(annotation_dict[i][4]) % 3 == 0:
        annotation_dict[i].append('inframe')
    else:
        annotation_dict[i].append('frameshift')

for i in sorted(annotation_dict):
    print annotation_dict[i]


header = 'Position' + '\t' + 'Ref' + '\t' + 'Alt' + '\t' + 'Type' + '\t' + 'Length' + '\t' + 'Indel (start/stop)' + \
         '\t' + 'Gene (start/stop)' + '\t' + 'Region (start/stop)' + '\t' + 'Position in Gene (start/stop)' + '\t' + \
         'Frame' + '\n'

with open(outfile, 'w') as outfile:
    outfile.write(header)
    to_write = []
    for i in sorted(annotation_dict):
        to_write.append('\t'.join(annotation_dict[i]))
    outfile.write('\n'.join(to_write))
