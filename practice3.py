from sys import argv
import os.path

script, infile1, outfile1 = argv

list1 = []
with open(infile1, 'r') as infile:
    for line in infile:
        if not line.startswith("Filename"):
            list1.append(line.strip().split("\t"))

count = 0

for item in list1:
        if int(item[1]) >= 1000 and item[-1].startswith("#"):
            count += 1

if os.path.isfile(outfile1):
    with open(outfile1, 'a') as outfile:
        outfile.write(str(infile1) + "\t")
        outfile.write(str(count)+"\n")
elif not os.path.isfile(outfile1):
    with open(outfile1, 'w') as outfile:
        outfile.write(str(infile1) + "\t")
        outfile.write(str(count)+"\n")
