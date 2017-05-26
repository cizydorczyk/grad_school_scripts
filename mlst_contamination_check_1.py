from sys import argv

script, input_vcf_samtoolsQ0, input_vcf_samtoolsQ13, output_full, output_g0, output_g5 = argv

class VcfPosition(object):
    def __init__(self, chromosome, position, ref, alt, quality, depth, dp4):
        self.chromosome = chromosome
        self.position = position
        self.ref = ref
        self.alt = alt
        self.quality = quality
        self.depth = depth
        self.dp4 = dp4

def parse_vcf(vcffile):
    vcf_objects = []
    with open(vcffile, 'r') as infile:
        for line in infile:
            if not line.startswith("#"):
                line1 = line.strip().split('\t')
                line1.pop(9)
                line1.pop(8)
                line1.pop(2)
                line1.pop(5)
                line2 = line1[-1].strip().split(";")
                dp = line2[0][3:]
                dp4 = None
                for i in line2:
                    if 'DP4' in i:
                        dp4 = line2[line2.index(i)][4:].split(",")
                vcf_objects.append(VcfPosition(line1[0],line1[1], line1[2], line1[3], line1[4], dp, dp4))
    return vcf_objects

Q0_vcf_objects_list = parse_vcf(input_vcf_samtoolsQ0)
Q13_vcf_objects_list = parse_vcf(input_vcf_samtoolsQ13)

def parse_variant_objects_list(objectslist):
    to_return = {}

    for i in objectslist:
        to_return[i.position] = [i.position]
        to_return[i.position] += [','.join(i.dp4)]
        ref_reads = float(i.dp4[0]) + float(i.dp4[1])
        alt_reads = float(i.dp4[2]) + float(i.dp4[3])
        percent_alt = str(round((alt_reads / (ref_reads + alt_reads) * 100), 2))
        to_return[i.position].append(percent_alt)
    return to_return

Q0_output = parse_variant_objects_list(Q0_vcf_objects_list)
Q13_output = parse_variant_objects_list(Q13_vcf_objects_list)

def parse_Q_output(Qoutputdict):
    outdict = {}
    for i in sorted(Qoutputdict):
        if Qoutputdict[i] > 0.0:
            outdict[i] = Qoutputdict[i]
    return outdict

Q0_final_dict = parse_Q_output(Q0_output)
Q13_final_dict = parse_Q_output(Q13_output)

for i in sorted(Q13_final_dict):
    Q13_final_dict[i].pop(0)
    Q0_final_dict[i] += Q13_final_dict[i]

for i in sorted(Q0_final_dict):
    if float(Q0_final_dict[i][2]) > 0.0 or float(Q0_final_dict[i][4] > 0.0):
        print Q0_final_dict[i]

header = 'Position' + '\t' + 'Q0_DP4' + '\t' + 'Q0_%_alt_reads' + '\t' + 'Q13_DP4' + 'Q13_%_alt_reads' + '\n'
with open(output_full, 'w') as outfile1:
    outfile1.write(header)
    to_write1 = []
    for i in sorted(Q0_final_dict):
        to_write1.append('\t'.join(Q0_final_dict[i]))
    outfile1.write('\n'.join(to_write1))

with open(output_g0, 'w') as outfile2:
    outfile2.write(header)
    to_write2 = []
    for i in sorted(Q0_final_dict):
        if float(Q0_final_dict[i][2]) > 0.0 or float(Q0_final_dict[i][4]) > 0.0:
            to_write2.append('\t'.join(Q0_final_dict[i]))
    outfile2.write('\n'.join(to_write2))

with open(output_g5, 'w') as outfile3:
    outfile3.write(header)
    to_write3 = []
    for i in sorted(Q0_final_dict):
        if float(Q0_final_dict[i][2]) >= 5.0 or float(Q0_final_dict[i][4]) >=5.0:
            to_write3.append('\t'.join(Q0_final_dict[i]))
    outfile3.write('\n'.join(to_write3))