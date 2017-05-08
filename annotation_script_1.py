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

def get_ranges(dict):
    gene_list = []
    for i in dict:
        gene_list.append(i)

    ranges = []
    for i in gene_list:
        range_ = set(range(*i))
        ranges.append(range_)
    return ranges

gene_ranges = get_ranges(TAB_ann_dict)
intergenic_ranges = get_ranges(GTF_intron_ann_dict)

for i in sorted(snps_dict):
    for range_ in gene_ranges:
        if i in range_:
            print i
            break
    else:
        for range_ in intergenic_ranges:
            if i in range_:
                print "intergenic"
                break


    # alt_base = snps_dict[i][1]
    # snp_pos = int(i)
    # sample_sequence = list(TAB_ann_dict[(483, 2027)].split('\t')[-2])
    # snp_index = (snp_pos - 483)
    #
    # codon = ''
    # if (snp_index + 1) % 3 == 0:
    #     codon = ''.join(sample_sequence[(snp_index-2):(snp_index + 1)])
    # elif (snp_index + 1) % 3 != 0 and (snp_index + 2) % 3 == 0:
    #     codon = ''.join(sample_sequence[(snp_index-1):(snp_index+2)])
    # elif (snp_index + 1) % 3 != 0 and (snp_index + 3) % 3 == 0:
    #     codon = ''.join(sample_sequence[snp_index:(snp_index + 3)])
    # print codon
    #
    # sample_sequence[snp_index] = alt_base
    # codon_alt = ''
    #
    # if (snp_index + 1) % 3 == 0:
    #     codon_alt = ''.join(sample_sequence[(snp_index-2):(snp_index + 1)])
    # elif (snp_index + 1) % 3 != 0 and (snp_index + 2) % 3 == 0:
    #     codon_alt = ''.join(sample_sequence[(snp_index-1):(snp_index+2)])
    # elif (snp_index + 1) % 3 != 0 and (snp_index + 3) % 3 == 0:
    #     codon_alt = ''.join(sample_sequence[snp_index:(snp_index + 3)])
    # print codon_alt
    #
    # codon_change = ''
    # if codons[codon] == codons[codon_alt]:
    #     codon_change = "synonymous"
    # elif codons[codon] != codons[codon_alt]:
    #     codon_change = "non-synonymous"
    # print codon_change
    # print i
