from sys import argv

script, tab_annotation, GTF_intron_annotation, snps_file = argv

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
            new_line = [line[2], line[3], line[4], str(int(line[4])-int(line[3]) + 1), line[8].split(';')[3], line[8].split(';')[4]]
            key = (int(line[3]), int(line[4]))
            GTF_intron_ann_dict[key] = '\t'.join(new_line)
# Genic annotation format:
# Sequence, Locus Tag, Feature Type, Start, End, Strand, Name, Product Name, Accession,
# GI, Length (nucleotides), MW (predicted), Length (amino acids), Nucleotide Sequence,
# Amino Acid Sequence
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
            TAB_ann_dict[key] = '\t'.join(line)

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

for i in sorted(snps_dict):
    bound = lookup_bounds(sorted_TAB_ann_bounds, i)
    if bound is not None:
        alt_base = snps_dict[i][1]
        snp_pos = int(i)
        sample_sequence = list(TAB_ann_dict[bound].split('\t')[-2])
        snp_index = (snp_pos - int(TAB_ann_dict[bound].split('\t')[3]))
        print '~~~~~'
        print bound
        print len(sample_sequence)
        print snp_index
        print '~~~~~~~'
        codon = ''
        if (snp_index + 1) % 3 == 0:
            codon = ''.join(sample_sequence[(snp_index-2):(snp_index + 1)])
        elif (snp_index + 1) % 3 != 0 and (snp_index + 2) % 3 == 0:
            codon = ''.join(sample_sequence[(snp_index-1):(snp_index+2)])
        elif (snp_index + 1) % 3 != 0 and (snp_index + 3) % 3 == 0:
            codon = ''.join(sample_sequence[snp_index:(snp_index + 3)])
        print codon

        sample_sequence[snp_index] = alt_base
        codon_alt = ''

        if (snp_index + 1) % 3 == 0:
            codon_alt = ''.join(sample_sequence[(snp_index-2):(snp_index + 1)])
        elif (snp_index + 1) % 3 != 0 and (snp_index + 2) % 3 == 0:
            codon_alt = ''.join(sample_sequence[(snp_index-1):(snp_index+2)])
        elif (snp_index + 1) % 3 != 0 and (snp_index + 3) % 3 == 0:
            codon_alt = ''.join(sample_sequence[snp_index:(snp_index + 3)])
        print codon_alt
        print i
        # codon_change = ''
        # if codons[codon] == codons[codon_alt]:
        #     codon_change = "synonymous"
        # elif codons[codon] != codons[codon_alt]:
        #     codon_change = "non-synonymous"
        # print codon_change
        # print i

    elif bound is None:
        print "intergenic"
