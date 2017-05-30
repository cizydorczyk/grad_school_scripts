from sys import argv

script, input_pileup = argv

class PileupRecord(object):
    def __init__(self, chromosome, position, reference, num_reads, reads, qualities, record):
        self.chromosome = chromosome
        self.position = position
        self.reference = reference
        self.num_reads = num_reads
        self.reads = reads
        self.qualities = qualities
        self.record = record

    def num_forward_ref_reads(self):
        reads = list(self.reads).count('.')
        return reads

    def num_reverse_ref_reads(self):
        reads = list(self.reads).count(',')
        return reads

    def num_alt_reads(self):
        reads_list = list(self.reads)
        # total_alt_reads = reads_list.count('C') + reads_list.count('c') + reads_list.count('G') + reads_list.count('g')\
        # + reads_list.count('A') + reads_list.count('a') + reads_list.count('T') + reads_list.count('t') + reads_list.count('N')\
        # + reads_list.count('n')
        bases_dict = {}
        bases_dict['C'] = reads_list.count('C') + reads_list.count('c')
        bases_dict['A'] = reads_list.count('A') + reads_list.count('a')
        bases_dict['T'] = reads_list.count('T') + reads_list.count('t')
        bases_dict['G'] = reads_list.count('G') + reads_list.count('g')
        return bases_dict



def parse_pileup(pileup_file):
    pileup_object_list = []
    with open(input_pileup, 'r') as infile:
        for line in infile:
            line = line.strip().split('\t')
            record = PileupRecord(line[0], line[1], line[2], float(line[3]), line[4], line[5], line)
            pileup_object_list.append(record)
    return pileup_object_list

pileup_objects = parse_pileup(input_pileup)

for i in pileup_objects:

    alt_dict = i.num_alt_reads()
    alt_sum = alt_dict['C'] + alt_dict['T'] + alt_dict['A'] + alt_dict['G']
    alt_counts = 'C:' + str(alt_dict['C']) + '\t' + 'A:' + str(alt_dict['A']) + '\t' + 'T:' + str(alt_dict['T']) + '\t' + 'G:' + str(alt_dict['G'])
    out_line = i.position + '\t' + str(i.num_reads) + '\t' + str(i.num_forward_ref_reads()) + '\t' + str(i.num_reverse_ref_reads())\
        + '\t' + alt_counts + '\t' + str(round(alt_sum / i.num_reads * 100,2))
    print out_line
