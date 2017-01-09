from sys import argv

# Deduplicated pos list does NOT have to exist prior to running this script.
script, totalposlist, dedupposlist = argv

dedup_list = []

infile = open(totalposlist, 'r')
outfile = open(dedupposlist, 'w')

temp3 = []
for line in infile:
    temp1 = line.split("\t")
    temp2 = []
    temp2.append(int(temp1[0]))
    temp2.append(str(temp1[1]))
    if temp2[0] not in dedup_list:
        temp3.append(temp2)
        dedup_list.append(temp2[0])

temp3.sort()

for i in temp3:
    outfile.write(str(i[0])+"\t")
    outfile.write(i[1])


infile.close()
outfile.close()

print "Hqposlist deduplicated and sorted."
