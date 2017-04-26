import itertools
import sys
import os.path

script, bwa, novo, last, bwa_list, novo_list, last_list = sys.argv

def parse_vcf(vcffile):
    positions = []
    with open(vcffile, 'r') as infile:
        for line in infile:
            if not line.startswith("#"):
                positions.append(line.strip().split('\t')[1])
    return positions


hq_bwa_positions = set(parse_vcf(bwa))
hq_novo_positions = set(parse_vcf(novo))
hq_last_positions = set(parse_vcf(last))

hq_positions = (hq_bwa_positions & hq_novo_positions & hq_last_positions)
print len(hq_positions)
