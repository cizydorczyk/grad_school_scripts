from sys import argv

script, list_of_isolate_snps, all_snps, output_1, output_2, output_3, output_4, output_5, output_6, output_7 = argv


isolate_snp_files = []
with open(list_of_isolate_snps, 'r') as infile1:
    for line in infile1:
        isolate_snp_files.append(line.strip())

all_snps_dict = {}

with open(all_snps, 'r') as infile2:
    for line in infile2:
        line_list = line.strip().split('\t')
        all_snps_dict[int(line_list[0])] = line_list[1]



##############################################
# # To count snps at ALL positions in reference genome:

all_pos_counts = {}
for i in range(1,6264405):
    all_pos_counts[i] = 0


for file_ in isolate_snp_files:
    isolate_snps = []
    with open(file_, 'r') as infile3:
        for line in infile3:
            if 'Position' not in line:
                isolate_snps.append(int(line.strip().split('\t')[0]))
    for i in isolate_snps:
        all_pos_counts[i] += 1

with open(output_1, 'w') as outfile1:
    to_write = []
    for i in sorted(all_pos_counts):
        to_write.append(str(i) + '\t' + str(all_pos_counts[i]))
    outfile1.write('\n'.join(to_write))

#############################################
# To count snps in bins of 1K:

tuple_list_1k = []
counts_dict_1k = {}

for i in range(0,6300):
    tuple_list_1k.append((i, i+1))
    counts_dict_1k[i+1] = 0

for file_ in isolate_snp_files:
    snps = []
    with open(file_, 'r') as infile3:
        for line in infile3:
            if 'Position' not in line:
                snps.append(int(line.strip().split('\t')[0]))

    print file_
    for i in snps:
        for min_, max_ in tuple_list_1k:
            if i > (max_*1000):
                continue
            elif (min_*1000) <= i:
                counts_dict_1k[max_] += 1
                break

print counts_dict_1k

with open(output_2, 'w') as outfile1:
    to_write = []
    for i in sorted(counts_dict_1k):
        to_write.append(str(i) + '\t' + str(counts_dict_1k[i]))
    outfile1.write('\n'.join(to_write))

##############################################
# To count snps in bins of 10K:

tuple_list_10k = []
counts_dict_10k = {}

for i in range(0,630):
    tuple_list_10k.append((i, i+1))
    counts_dict_10k[i+1] = 0

for file_ in isolate_snp_files:
    snps = []
    with open(file_, 'r') as infile3:
        for line in infile3:
            if 'Position' not in line:
                snps.append(int(line.strip().split('\t')[0]))

    print file_
    for i in snps:
        for min_, max_ in tuple_list_10k:
            if i > (max_*10000):
                continue
            elif (min_*10000) <= i:
                counts_dict_10k[max_] += 1
                break

print counts_dict_10k

with open(output_3, 'w') as outfile1:
    to_write = []
    for i in sorted(counts_dict_10k):
        to_write.append(str(i) + '\t' + str(counts_dict_10k[i]))
    outfile1.write('\n'.join(to_write))

##############################################
# To count snps in bins of 100K:

tuple_list_100k = []
counts_dict_100k = {}

for i in range(0,63):
    tuple_list_100k.append((i, i+1))
    counts_dict_100k[i+1] = 0

for file_ in isolate_snp_files:
    snps = []
    with open(file_, 'r') as infile3:
        for line in infile3:
            if 'Position' not in line:
                snps.append(int(line.strip().split('\t')[0]))

    print file_
    for i in snps:
        for min_, max_ in tuple_list_100k:
            if i > (max_*100000):
                continue
            elif (min_*100000) <= i:
                counts_dict_100k[max_] += 1
                break

print counts_dict_100k

with open(output_4, 'w') as outfile1:
    to_write = []
    for i in sorted(counts_dict_100k):
        to_write.append(str(i) + '\t' + str(counts_dict_100k[i]))
    outfile1.write('\n'.join(to_write))

##############################################
# To count distribution of snps in reference genome:
# Note this is not the same as below; this gets the number mutated sites per
# X K of the reference genome, but only counts each mutated position ONCE,
# rather than getting the total number of isolates with a SNP at each site.

tuple_list1k = []
counts_dict_1k = {}

for i in range(0,6300):
    tuple_list1k.append((i, i+1))
    counts_dict_1k[i+1] = 0

tuple_list10k = []
counts_dict_10k = {}

for i in range(0,630):
    tuple_list10k.append((i, i+1))
    counts_dict_10k[i+1] = 0

tuple_list100k = []
counts_dict_100k = {}

for i in range(0,63):
    tuple_list100k.append((i, i+1))
    counts_dict_100k[i+1] = 0

for i in snps:
    for min_, max_ in tuple_list_1k:
        if i > (max_*1000):
            continue
        elif (min_*1000) <= i:
            counts_dict_1k[max_] += 1
            break

for i in snps:
    for min_, max_ in tuple_list_10k:
        if i > (max_*10000):
            continue
        elif (min_*10000) <= i:
            counts_dict_10k[max_] += 1
            break

for i in snps:
    for min_, max_ in tuple_list_100k:
        if i > (max_*100000):
            continue
        elif (min_*100000) <= i:
            counts_dict_100k[max_] += 1
            break

with open(output_5, 'w') as outfile1:
    to_write = []
    for i in sorted(counts_dict_1k):
        to_write.append(str(i) + '\t' + str(counts_dict_1k[i]))
    outfile1.write('\n'.join(to_write))

with open(output_6, 'w') as outfile1:
    to_write = []
    for i in sorted(counts_dict_10k):
        to_write.append(str(i) + '\t' + str(counts_dict_10k[i]))
    outfile1.write('\n'.join(to_write))

with open(output_7, 'w') as outfile1:
    to_write = []
    for i in sorted(counts_dict_100k):
        to_write.append(str(i) + '\t' + str(counts_dict_100k[i]))
    outfile1.write('\n'.join(to_write))
