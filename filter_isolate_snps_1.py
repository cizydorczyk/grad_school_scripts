from sys import argv
import os.path

script, input_snp_file, removed_positions_file, common_log_file, snp_output_file = argv

isolate_number = input_snp_file.split('/')[-1].split('_')[0]
print isolate_number

snp_dict = {}
with open(input_snp_file, 'r') as infile1:
    for line in infile1:
        if not line.startswith("Position"):
            position = int(line.strip().split('\t')[0])
            if position not in snp_dict:
                snp_dict[position] = [line.strip().split('\t')[1], line.strip().split('\t')[2]]

to_remove = []
with open(removed_positions_file, 'r') as infile2:
    for line in infile2:
        to_remove.append(int(line.strip()))

original_length = len(snp_dict)
print "original length", original_length
print "total num positions to remove ", len(to_remove)

count = 0
# For each position to remove, if it appears in the isolate's snps, remove that snp:
for i in to_remove:
    if i in snp_dict:
        count += 1
        del snp_dict[i]

print "num snps removed ", count

ref_calls = []
for i in snp_dict:
    if snp_dict[i][0] == snp_dict[i][1]:
        ref_calls.append(i)

count2 = 0
for i in ref_calls:
    count2 += 1
    del snp_dict[i]

print "num ref calls removed", count2

print "new length ", len(snp_dict)

if os.path.isfile(common_log_file):
    with open(common_log_file, 'a') as outfile1:
        outfile1.write("Isolate " + isolate_number + '\n' +
                       "original length " + str(original_length) + '\n' +
                       "total possible num snp positions to remove " + str(len(to_remove)) + '\n' +
                       "num snps removed " + str(count) + '\n' +
                       "num reference calls removed " + str(count2) + '\n' +
                       "new length " + str(len(snp_dict)) + '\n' + '\n')

elif not os.path.isfile(common_log_file):
    with open(common_log_file, 'w') as outfile1:
        outfile1.write("Isolate " + isolate_number + '\n' +
                       "original length " + str(original_length) + '\n' +
                       "total possible num snp positions to remove " + str(len(to_remove)) + '\n' +
                       "num snps removed " + str(count) + '\n' +
                       "num reference calls removed " + str(count2) + '\n' +
                       "new length " + str(len(snp_dict)) + '\n' + '\n')

with open(snp_output_file, 'w') as outfile2:
    outfile2.write("Position" + '\t' + "Reference" + '\t' + "Base_called" + '\n')
    for i in sorted(snp_dict):
        outfile2.write(str(i) + '\t' + '\t'.join(snp_dict[i]) + '\n')