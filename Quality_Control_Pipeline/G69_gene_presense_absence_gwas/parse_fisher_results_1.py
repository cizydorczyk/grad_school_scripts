from sys import argv

script, fisherinput = argv

class FisherOutput(object):
    def __init__(self, snp, oddsratio, pvalue):
        self.snp = snp
        self.oddsratio = oddsratio
        self.pvalue = pvalue

fisher_objects = []
undef_count = 0
def_count = 0
undef_snps = []
def_snps = []
with open(fisherinput, 'r') as infile1:
    for line in infile1:
        temp_line = line.strip().split('\t')
        fisher_objects.append(FisherOutput(temp_line[0], temp_line[1], temp_line[2]))

for i in fisher_objects:
    if i.oddsratio == '0.0':
        undef_count += 1
        undef_snps.append(i)
    elif i.oddsratio == 'inf':
        undef_count += 1
        undef_snps.append(i)
    else:
        def_count += 1
        def_snps.append(i)

print undef_count, def_count

sig_count = 0
unsig_count = 0
bf_pvalue = 0.05 / float(def_count)

for i in def_snps:
    # print i.pvalue, bf_pvalue

    if i.pvalue <= bf_pvalue:
        sig_count += 1
    else:
        unsig_count += 1

print sig_count, unsig_count
