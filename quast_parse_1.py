from sys import argv
import os.path

script, isolate_number, inputfile1, outputtsv = argv

with open(inputfile1, 'r') as infile:
    list1 = []
    for line in infile:
        list1.append(line.strip('\n').split('\t'))

def better_assembly(isolatenumber, tsvlist):
    qc5_count = 0
    qc20_count = 0
# Total length (>= 0 bp)
    if tsvlist[7][1] > tsvlist[7][2]:
        qc5_count += 1
    elif tsvlist[7][1] < tsvlist[7][2]:
        qc20_count += 1
    else:
        pass
# Total length (>= 1000 bp)
    if tsvlist[8][1] > tsvlist[8][2]:
        qc5_count += 1
    elif tsvlist[8][1] < tsvlist[8][2]:
        qc20_count += 1
    else:
        pass
# Total length (>= 5000 bp)
    if tsvlist[9][1] > tsvlist[9][2]:
        qc5_count += 1
    elif tsvlist[9][1] < tsvlist[9][2]:
        qc20_count += 1
    else:
        pass
# Total length (>= 10,000 bp)
    if tsvlist[10][1] > tsvlist[10][2]:
        qc5_count += 1
    elif tsvlist[10][1] < tsvlist[10][2]:
        qc20_count += 1
    else:
        pass
# Total length (>= 25,000 bp)
    if tsvlist[11][1] > tsvlist[11][2]:
        qc5_count += 1
    elif tsvlist[11][1] < tsvlist[11][2]:
        qc20_count += 1
    else:
        pass
# Total length (>= 50,000 bp)
    if tsvlist[12][1] > tsvlist[12][2]:
        qc5_count += 1
    elif tsvlist[12][1] < tsvlist[12][2]:
        qc20_count += 1
    else:
        pass
# Number of contigs
    if tsvlist[13][1] > tsvlist[13][2]:
        qc20_count += 1
    elif tsvlist[13][1] < tsvlist[13][2]:
        qc5_count += 1
    else:
        pass
# Largest contig
    if tsvlist[14][1] > tsvlist[14][2]:
        qc5_count += 1
    elif tsvlist[14][1] < tsvlist[14][2]:
        qc20_count += 1
    else:
        pass
# Total length (in contigs > 500 bp)
    if tsvlist[15][1] > tsvlist[15][2]:
        qc5_count += 1
    elif tsvlist[15][1] < tsvlist[15][2]:
        qc20_count += 1
    else:
        pass
# N50
    if tsvlist[17][1] > tsvlist[17][2]:
        qc5_count += 1
    elif tsvlist[17][1] < tsvlist[17][2]:
        qc20_count += 1
    else:
        pass
# N75
    if tsvlist[18][1] > tsvlist[18][2]:
        qc5_count += 1
    elif tsvlist[18][1] < tsvlist[18][2]:
        qc20_count += 1
    else:
        pass
# L50
    if tsvlist[19][1] > tsvlist[19][2]:
        qc20_count += 1
    elif tsvlist[19][1] < tsvlist[19][2]:
        qc5_count += 1
    else:
        pass
# L75
    if tsvlist[20][1] > tsvlist[20][2]:
        qc20_count += 1
    elif tsvlist[20][1] < tsvlist[20][2]:
        qc5_count += 1
    else:
        pass
# N's per 100 kbp
    if tsvlist[21][1] > tsvlist[21][2]:
        qc20_count += 1
    elif tsvlist[21][1] < tsvlist[21][2]:
        qc5_count += 1
    else:
        pass
    result = str(isolatenumber) + '\t' + str(qc5_count) + "\t" + str(qc20_count)
    return result

result2 = better_assembly(isolate_number, list1)
print "Now parsing: " + inputfile1
print result2

if os.path.isfile(outputtsv):
    with open(outputtsv, 'a') as outfile:
        outfile.write(result2 + '\n')
elif not os.path.isfile(outputtsv):
    with open(outputtsv, 'w') as outfile:
        outfile.write("isolate_number" + '\t' + '# best in qc5' + '\t' + '# best in qc 20' + '\n')
        outfile.write(result2 + '\n')
