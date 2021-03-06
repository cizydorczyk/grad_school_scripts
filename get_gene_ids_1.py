from sys import argv
import os.path

script, input_annotation, output_gene_list, isolate_records_output = argv

class GenicAnnotationObject(object):
    def __init__(self, snp_pos, ref_base, alt_base, pos_in_gene, ref_codon, alt_codon, pos_in_codon, ref_aa, alt_aa, \
                 type_, sequence, locus_tag, feature_type, start, end, strand, gene_name, product_name, accession, GI, \
                 length_nuc, MW, length_aa):
        self.snp_pos = snp_pos
        self.ref_base = ref_base
        self.alt_base = alt_base
        self.pos_in_gene = pos_in_gene
        self.ref_codon = ref_codon
        self.alt_codon = alt_codon
        self.pos_in_codon = pos_in_codon
        self.ref_aa = ref_aa
        self.alt_aa = alt_aa
        self.type = type_
        self.sequence = sequence
        self.locus_tag = locus_tag
        self.feature_type = feature_type
        self.start = start
        self.end = end
        self.strand = strand
        self.gene_name = gene_name
        self.product_name = product_name
        self.accession = accession
        self.GI = GI
        self.length_nuc = length_nuc
        self.MW = MW
        self.length_aa = length_aa

gene_annotations = []
with open(input_annotation, 'r') as infile:
    for line in infile:
        if not line.startswith("#") and 'intergenic' not in line:
            annotation = line.strip().split("\t")
            gene_annotations.append(GenicAnnotationObject(annotation[0], annotation[1], annotation[2], annotation[3], \
                                                          annotation[4], annotation[5], annotation[6], annotation[7], \
                                                          annotation[8], annotation[9], annotation[10], annotation[11],\
                                                          annotation[12], annotation[13], annotation[14], annotation[15]\
                                                          , annotation[16], annotation[17], annotation[18], annotation[19]\
                                                          , annotation[20], annotation[21], annotation[22]))

print 'number of genic annotations', len(gene_annotations)

known_gene_ids = []
unknown_gene_ids = []

for i in gene_annotations:
    if 'NA' not in i.accession:
        known_gene_ids.append(i.accession)
    else:
        unknown_gene_ids.append(i.accession)

known_gene_ids_set = set(known_gene_ids)

isolate_record = str(input_annotation.split('/')[-1].split('_')[0]) + '\t' + str(len(gene_annotations)) + '\t' + \
    str(len(known_gene_ids)) + '\t' + str(len(unknown_gene_ids)) + '\t' + str(len(known_gene_ids_set))
print isolate_record

with open(output_gene_list, 'w') as outfile:
    outfile.write('\n'.join(known_gene_ids_set))

if not os.path.isfile(isolate_records_output):
    with open(isolate_records_output, 'w') as outfile2:
        outfile2.write('Isolate' + '\t' + 'Tot # Genic SNPs' + '\t' + '# Known Acc' + '\t' +
                       '# Unknown Acc' + '\t' + '# Unique Gene Acc' + '\n')
        outfile2.write(isolate_record + '\n')

elif os.path.isfile(isolate_records_output):
    with open(isolate_records_output, 'a') as outfile2:
        outfile2.write(isolate_record + '\n')