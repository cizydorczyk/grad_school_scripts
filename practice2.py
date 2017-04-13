from sys import argv

script, inputfile = argv


list1 = []
with open(inputfile, 'r') as infile:
    for line in infile:
        list1.append(line.strip().split("\t"))
list1[0].pop(8)
print '\t'.join(list1[0])
for line in list1:
    if "ST" not in line:
        if len(line) == 8:
            print '\t'.join(line)
        else:
            line.pop(8)
            print '\t'.join(line)
