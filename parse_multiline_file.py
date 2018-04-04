# Eg. for multiline fasta or fastq...also includes method for when end of file
# is reached (except IndexError part)

from sys import argv

script, inputfile = argv

with open(inputfile, 'r') as infile:
    list1 = list(infile)

dict1 = {}
for num, i in enumerate(list1):
    if i.startswith(">"):
        #print i,
        tick = 1
        sequence = ''
        j = list1[num+tick]
        while not j.startswith(">"):
            sequence += j.strip()
            tick += 1
            try:
                j = list1[num+tick]
            except IndexError:
                break
        dict1[i] = sequence

print dict1
