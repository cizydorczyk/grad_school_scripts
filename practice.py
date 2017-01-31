from sys import argv
from oset import oset
script, inputfile = argv

list1 = []
with open(inputfile, 'r') as infile:
    for line in infile:
        if "orange" in line:
            line1 = line.split("\t")
            list1.append(line1[0])

list2 = []
for line in list1:
    line2 = line.split("_")
    list2.append(int(line2[-1]))
list3 = range(1,900)
print len(list3)
print len(list2)
set2 = oset(list2)
set3 = oset(list3)

count = 0
for i in set3:
    if i not in set2:
        print i,
        count += 1
print "Total count is :" + str(count)
